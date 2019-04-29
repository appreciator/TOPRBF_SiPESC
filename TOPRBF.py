# -*- coding: utf-8 -*-
'''
===============================================================================
TOPRBF
Version 1.0, December 16, 2016
Python program developed for the 2D mean compliance minimum problem
with the parameterized level set method using radial basis functions.
Developed by: Zuyu Li, Peng Wei
School of Civil Engineering and Transportation
South China University of Technology 
Email: lizuyu0091@163.com
===============================================================================
'''
from __future__ import division #整数相除的浮点数运算
from numpy import *
from scipy import sparse
from scipy.sparse.linalg import spsolve
import datetime
import time
from PyQt5.QtCore import QDir
from PyQt5.QtCore import QTime
from matplotlib import pyplot
pyplot.ion() #matplotlib绘图交互模式
startTime = datetime.datetime.now() #开始计时

## SiPESC ##
#模型文件及路径设定
workSpace = 'D:\\SiPESC_WorkSpace\TestModel' #工作空间路径
programName = 'DB' #数据库文件夹名称
modelName = 'Model_'+str(QTime().currentTime().toString()) #模型数据文件夹名称
modelName = modelName.replace(":","_")
print 'Model: ',modelName
fileName = 'CantileverBeam.cdb' #导入模型文件名

objectManager = MObjectManager()
extManager = MExtensionManager()
dbManager = objectManager.getObject('org.sipesc.core.engdbs.mdatabasemanager')
vecFactory = objectManager.getObject('org.sipesc.fems.matrix.vectorfactory') #MVector
matFactory = objectManager.getObject('org.sipesc.fems.matrix.matrixfactory') #MMatrix
controlTrans = extManager.createExtension('org.sipesc.fems.controlmatrix.MControlMatrixTrans')
mcTrans = extManager.createExtension('org.sipesc.fems.scriptassist.MMatrixAssistor') #含multiply函数用于MMatrix与MVector或标量的相乘,可转换单刚类型为MMatrix

#数据库初始化
QDir.setCurrent(workSpace)
option = dbManager.createDatabaseOption()
option.setAddress(workSpace)
option.setCacheName('org.sipesc.core.engdbs.data.lrucache')
option.setDeviceName('org.sipesc.core.engdbs.data.file')
option.setDatabaseName(programName)
db = dbManager.createDatabase('org.sipesc.engdbs') #创建数据库
isOK = db.open(option,False) #根据设定的选项打开数据库 True:只读 False:可读可写
db.clear() #清空初始化的数据库

#模型数据初始化
model = dbManager.createDataModel()
isOK = model.open(db,modelName) #在数据库中创建模型数据空间
model.clear() #清空初始化的模型数据文件夹

#从第三方文件中导入模型数据到数据库
factoryManager = extManager.createExtension('org.sipesc.utilities.MExtensionFactoryManager')
factoryManager.initialize('fems.factory.dataimportor')
factory = factoryManager.getFactory('cdb') #支持读取: bdf, cdb, inp
importor = factory.createExtension()
monitor = importor.getProgressMonitor()
indicator = dbManager.createIndicator(monitor) 
importor.initialize(model)
importor.import_(fileName) #模型文件导入

#定义优化循环外的前处理流程--含原始单刚计算
pre_static_Commands=['org.sipesc.fems.femstask.MNodeMapManager',
'org.sipesc.fems.femstask.MNodeSortBySpectrum',
'org.sipesc.fems.controlmatrix.MDofStandardParserManager',
'org.sipesc.fems.controlmatrix.MDofNumberingStandardManager',
'org.sipesc.fems.controlmatrix.MNodeControlMatrixStandardManager',
'org.sipesc.fems.controlmatrix.MElementLocalCoordTransStandardManager',
'org.sipesc.fems.controlmatrix.MElementControlMatrixStandardManager',
'org.sipesc.fems.femstask.MSolidElementStiffManager']

#定义优化循环内的荷载与反力计算(依赖于更改后的单刚)、总刚组装、方程组求解等流程
solve_static_Commands=['org.sipesc.fems.femstask.MGivenValueLoadVectorManager',
'org.sipesc.fems.femstask.MLoadComponentManager',
'org.sipesc.fems.femstask.MLoadManager',
'org.sipesc.fems.femstask.MResultsManager']

#定义优化循环外的后处理流程
outpost_static_Commands = ['org.sipesc.fems.strain.MSolidElementStrainManager',
'org.sipesc.fems.stress.MSolidElementStressManager',
'org.sipesc.fems.stress.MSolidNodeStressManager']

def remove_cfg(model):
    #清除数据库中各类文件
    model.remove("GlobalStiffDiag",2)
    model.remove("GlobalStiffMatrix",4)


#有限元前处理各流程依次计算
localCoords = dbManager.createDataManager()
isOK = localCoords.open(model,'LocalCoordTransMatrix',False)
localCoords.clear() #清空局部坐标数据
indicator = extManager.createExtension('org.sipesc.utilities.MSimpleProgressIndicator')
commands = pre_static_Commands
nCommand = len(commands)
for iCom in range(1,nCommand+1,1):
	task = extManager.createExtension(commands[iCom-1])
	task.initialize(model)
	monitor = task.getProgressMonitor()
	indicator.setMonitor(monitor)
	task.start()

#从数据库中获取有限元模型信息
eleControlManager = dbManager.createDataManager()
isOK = eleControlManager.open(model,'ElementControlMatrix',True)
eleKManager = dbManager.createDataManager()
isOK = eleKManager.open(model,'ElementStiffMatrix',False) #在模型数据文件夹中打开单刚文件
nEle = eleKManager.getDataCount() #单元数目
#构建单元的节点组成矩阵--所用SiPESC版本不支持打开model里的model,故暂时人工构造
#矩形结构,边长为1的正方形单元,左下角节点号为1,坐标(0,0),右上角节点号为nNode,坐标(nelx,nely),按Y正向优先排序
nelx = 60 #X方向单元数
nely = 30 #Y方向单元数
eleN1 = tile(arange(1,nely+1,1),(nelx,1)).T+kron(arange(0,nelx,1),(nely+1)*ones((nely,1),dtype=int32)) #节点号从1开始,整数变量
elementNode = tile(eleN1.ravel('F').reshape(eleN1.size,1),(1,4))+tile(array([0,nely+1,nely+2,1]),(nelx*nely,1))
'''
elePath = dbManager.createDataModel()
isOK = elePath.open(model,'ElementPath')
elementNode = zeros((nEle,4)) #仅考虑全为四节点单元的情况
for iEle in range(1,nEle+1,1):
	ieleData = elePath.find(iEle)
	eleDataManager = ieleData.getCurrentManager()
	eleData = eleDataManager.getData(iEle)
	for ieleNode in range(eleData.NodeCount)
		elementNode[iEle-1,ieleNode] = eleData.getNodeId(ieleNode)	
'''
#构建节点坐标矩阵
nodeCoorManager = dbManager.createDataManager()
isOK = nodeCoorManager.open(model,'Node',True)
nNode = nodeCoorManager.getDataCount() #节点数目
nodeCoor = zeros((nNode,2))
for iNode in range(1,nNode+1,1):
	nodeData = nodeCoorManager.getDataAt(iNode-1)
	nodeCoor[iNode-1,:] = [nodeData.getX(),nodeData.getY()]
#各单元的初始材料参数
materialManager = dbManager.createDataManager()
isOK = materialManager.open(model,'Material',False)
nMat = materialManager.getDataCount() #材料类型数目
#Plane42-Plane1042:密度,弹模,泊松比
eleDensity = zeros((nEle,1)) #各单元材料密度
eleYoungE = zeros((nEle,1)) #各单元弹模
if (nMat==nEle):
	#考虑自重时,材料类型数应等于单元数
	for iMat in range(nMat):
		eleMat = materialManager.getDataAt(iMat)
		eleDensity[iMat,0] = eleMat.getValue(0)
		eleYoungE[iMat,0] = eleMat.getValue(1)	
else:
	for iEle in range(nEle):
		#暂无法打开ElementPath,故此段代码未测试
		'''
		ieleData = elePath.find(iEle)
		eleDataManager = ieleData.getCurrentManager()
		eleData = eleDataManager.getData(iEle)	
		iMat = eleData.PropertyId
		eleMat = materialManager.getDataAt(iMat)
		'''
		eleMat = materialManager.getDataAt(0) #暂假定各单元均为第一种材料
		eleDensity[iEle,0] = eleMat.getValue(0)
		eleYoungE[iEle,0] = eleMat.getValue(1)	
	
## 水平集函数初始化 ##
Phi = zeros((nNode,1)) #初始结构为全设计域实体,亦可选用以下针对meshgrid型数据的孔洞初始化方法
'''
r = nely*0.1 #孔洞半径
hX = nelx*column_stack([tile(array([1.0/6,5.0/6]),(1,3)),tile(array([0,1.0/3,2.0/3,1]),(1,2)),1.0/2])
hY = nely*column_stack([kron(array([0,1.0/2,1]),ones((1,2))),kron(array([1.0/4,3.0/4]),ones((1,4))),1.0/2])
X,Y = meshgrid(linspace(0,nelx,nelx+1),linspace(0,nely,nely+1))
dX = tile(X,(hX.size,1,1))-reshape(hX,(hX.size,1,1))
dY = tile(Y,(hY.size,1,1))-reshape(hY,(hY.size,1,1))
Phi = sqrt(dX**2+dY**2)-r
Phi = maximum(-3,minimum(3,Phi.min(0))).reshape((-1,1),order='F')
'''
## 径向基函数参数化 ##
cRBF = 1e-4 #RBF参数
Ax = nodeCoor[:,0].reshape(-1,1)- nodeCoor[:,0]
Ay = nodeCoor[:,1].reshape(-1,1)- nodeCoor[:,1]
A = sqrt(Ax**2+Ay**2+cRBF**2)
G = hstack([hstack([A,ones((nNode,1))]),hstack([nodeCoor[:,0].reshape(-1,1),nodeCoor[:,1].reshape(-1,1)])])
G = vstack([G,hstack([vstack([ones((1,nNode)),vstack([nodeCoor[:,0],nodeCoor[:,1]])]),zeros((3,3))])])
pGpX = vstack([hstack([1.0*Ax/A,tile(array([0,1,0]),(nNode,1))]),hstack([tile(array([0,1,0]),(nNode,1)).T,zeros((3,3))])])
pGpY = vstack([hstack([1.0*Ay/A,tile(array([0,0,1]),(nNode,1))]),hstack([tile(array([0,0,1]),(nNode,1)).T,zeros((3,3))])])
Alpha = linalg.solve(G,row_stack([Phi,0,0,0]))
## 优化循环  ##
volfrac = 0.5 #体积分数
nLoop = 200 #最大循环步数
nRelax = 30 #体积松弛步数
dt = 1 #时间演化步长
delta = 10
mu = 0.1
gama = 0.05
comp = zeros((nLoop,1))
vol = zeros((nLoop,1))
eleComp = zeros((nEle,1))
equivE = zeros((nEle,nLoop+1))
equivE[:,0] = eleYoungE.copy().ravel('F') #原始全域实体模型的弹模
for iT in range(1,nLoop+1,1):   
	print 'Iteration number: ',iT 	
	#清空部分数据
	remove_cfg(model)
	#结果向量数据清空	
	resultManager = dbManager.createDataManager()
	isOK = resultManager.open(model,'ResultVector',False)	
	resultManager.clear()	
	#各单元实体部分体积分数计算
	s,t = meshgrid(linspace(-1,1,21),linspace(-1,1,21)) 
	tmpPhi1 = dot(((1-s.ravel('F'))*(1-t.ravel('F'))/4).reshape(-1,1),Phi[elementNode[:,0]-1,0].reshape(1,-1))
	tmpPhi2 = dot(((1+s.ravel('F'))*(1-t.ravel('F'))/4).reshape(-1,1),Phi[elementNode[:,1]-1,0].reshape(1,-1))
	tmpPhi3 = dot(((1+s.ravel('F'))*(1+t.ravel('F'))/4).reshape(-1,1),Phi[elementNode[:,2]-1,0].reshape(1,-1))
	tmpPhi4 = dot(((1-s.ravel('F'))*(1+t.ravel('F'))/4).reshape(-1,1),Phi[elementNode[:,3]-1,0].reshape(1,-1))
	tmpPhi = tmpPhi1+tmpPhi2+tmpPhi3+tmpPhi4
	eleVol = 1.0*(tmpPhi>=0).sum(axis=0)/s.size #单元实体材料体积分数
	equivE[:,iT] = (eleYoungE.ravel('F')*1e-9+eleVol*(eleYoungE.ravel('F')-eleYoungE.ravel('F')*1e-9))
	vol[iT-1,0] = eleVol.mean() #结构总体积分数, 此处假定各单元体积相同
	print 'Volume: ',vol[iT-1,0]	
	#修改数据库中的材料
	if (nMat==nEle):
		#仅当考虑自重,即材料数等于单元数时,对材料数据进行修改
		for iMat in range(nMat):	
			eleMat = materialManager.getDataAt(iMat)
			eleMat.setValue(0,eleVol[iMat]*eleDensity[iMat,0]) #修改密度
			eleMat.setValue(1,equivE[iMat,iT]) #修改弹模	
	
	#修改数据库中的单刚
	for iEle in range(1,nEle+1,1):	
		eleKData = eleKManager.getData(iEle) #提取当前单元的单刚
		eleK = mcTrans.convertStiffMatrix(eleKData)
		tmp_eleK = mcTrans.multiply(eleK,equivE[iEle-1,iT]/equivE[iEle-1,iT-1]) #按等效弹模缩放单刚			
		tmp_eleKData = tmp_eleK.toStiffData()	
		tmp_eleKData.setId(iEle)
		eleKManager.setData(tmp_eleKData)	

	#有限元荷载计算、总刚组装、方程组求解
	commands = solve_static_Commands
	nCommand = len(commands)
	for iCom in range(1,nCommand+1,1):
		task = extManager.createExtension(commands[iCom-1])
		task.initialize(model)
		monitor = task.getProgressMonitor()
		indicator.setMonitor(monitor)
		task.start()	
		
	#应变能计算
	resultManager = dbManager.createDataManager()
	isOK = resultManager.open(model,'ResultVector',True)
	dispData = resultManager.getDataAt(0)
	dispVec = vecFactory.createVector()
	dispVec.fromData(dispData) #提取位移向量	
	calEleUKU = extManager.createExtension('org.sipesc.fems.femstask.MComputeElementKU')
	for iEle in range(1,nEle+1,1):	
		eleKData = eleKManager.getData(iEle) #提取当前单元的单刚	
		eleControlData = eleControlManager.getData(iEle) #提取当前单元的单元控制矩阵				
		eleComp[iEle-1,0] = calEleUKU.computeElementUTKU(eleControlData,eleKData,dispVec)

	comp[iT-1,0] = eleComp.sum()
	print 'Compliance: ',comp[iT-1,0]	

	## matplotlib可视化 ##
	#收敛曲线
	pyplot.figure(1,facecolor='white')
	pyplot.subplot(2,1,1)
	pyplot.plot(arange(1,iT+1,1),comp[0:iT,0],'r*-')	
	pyplot.xlabel('Iteration number')
	pyplot.ylabel('Compliance')
	pyplot.subplot(2,1,2)
	pyplot.plot(arange(1,iT+1,1),vol[0:iT,0],'bo-')
	pyplot.xlabel('Iteration number')
	pyplot.ylabel('Volume fraction')
	pyplot.show()
	time.sleep(1)
	pyplot.close()
	#拓扑结果--仅适用于meshgrid型数据
	fig = pyplot.figure(2,facecolor='white')
	ax = fig.add_subplot(1,1,1)
	ax.axis('equal')
	ax.axis('scaled')	
	ax.set_xticks([])
	ax.set_yticks([])		
	X = nodeCoor[:,0].reshape((nely+1,nelx+1),order='F')
	Y = nodeCoor[:,1].reshape((nely+1,nelx+1),order='F')
	CS = ax.contourf(X,Y,Phi.reshape((nely+1,nelx+1),order='F'),[0,0],cmap=pyplot.cm.gray,zorder=1)
	pyplot.xlim(nodeCoor[:,0].min(),nodeCoor[:,0].max())
	pyplot.ylim(nodeCoor[:,1].min(),nodeCoor[:,1].max())
	CS2 = pyplot.contour(CS,levels=CS.levels,linestyles='solid',linewidths=2,colors='r',zorder=2) 
	pyplot.show()	
	time.sleep(2)
	pyplot.close()
	## 收敛判断 ##
	if (iT>nRelax) & (abs(vol[iT-1,0]-volfrac)/volfrac<1e-3) & (
	all(abs(comp[iT-1,0]-comp[iT-10:iT-2,0])/comp[iT-1,0]<1e-3)):		
		break
		
	## 拉格朗日乘子 ##	
	if (iT <= nRelax):
		lag = 20*mu*maximum((vol[iT-1,0]-vol[0,0]+(vol[0,0]-volfrac)*iT/nRelax),0)
	else:
		lag += gama*(vol[iT-1,0]-volfrac)
		gama = minimum(gama+0.05,1)	  

	## 水平集函数演化 ##
	gradPhi = sqrt(dot(pGpX,Alpha)**2+dot(pGpY,Alpha)**2)	
	indexDelta = abs(Phi)<=delta
	DeltaPhi = zeros(Phi.shape)	
	DeltaPhi[indexDelta] = 0.75/delta*(1-Phi[indexDelta]**2/delta**2)
	nodeComp = zeros((nNode,1))
	nodeCompCoeff = zeros((nNode,1))
	for iEle in range(nEle):
		eleNode = elementNode[iEle,:]
		nodeComp[eleNode-1,0] += eleComp[iEle,0]/1 #假定各单元体积相同
		nodeCompCoeff[eleNode-1,0] += 1
	nodeComp /= nodeCompCoeff #节点应变能密度	
	B = (nodeComp.ravel('F')/nodeComp.mean()-lag)*DeltaPhi.ravel('F')*delta/0.75	
	Alpha += dt*linalg.solve(G,row_stack([B.reshape(-1,1),0,0,0]))	
	if any((eleVol<1) & (eleVol>0)):
		Alpha /= mean(gradPhi[unique(elementNode[((eleVol<1) & (eleVol>0)),:])])
	Phi = dot(G[0:-3,:],Alpha)			

#循环外的有限元后处理
commands = outpost_static_Commands
nCommand = len(commands)
for iCom in range(1,nCommand+1,1):
	task = extManager.createExtension(commands[iCom-1])
	task.initialize(model)
	monitor = task.getProgressMonitor()
	indicator.setMonitor(monitor)
	task.start()	
	
endTime = datetime.datetime.now()
runTime = (endTime-startTime).seconds
print 'Elapsed time(s): ',runTime

'''
===============================================================================
Reference:
Zuyu Li, Peng Wei.(2016).An 88-line MATLAB code for the parameterized level set
method based topology optimization using radial basis functions

***********************************Disclaimer***********************************
The authors reserve all rights for the program. The programs may be distributed 
and used for academic and educational purposes. The authors do not guarantee that
the code is free from errors, and they shall not be liable in any event caused by
the use of the program.
===============================================================================
'''