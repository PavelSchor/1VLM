import vtk
import numpy as np

def uniqueRows(a):
	a = np.ascontiguousarray(a)
	unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
	return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))


def PGtoVTK(fname,pg,plotWake=False):
	
	def pgToPolyData(pg):
		points = vtk.vtkPoints()
		cells = vtk.vtkCellArray()


		arr_cp=vtk.vtkFloatArray()
		arr_cp.SetName('CP')
		arr_pload=vtk.vtkFloatArray()
		arr_pload.SetName('PLOAD')
		arr_eids=vtk.vtkFloatArray()
		arr_eids.SetName('GID')
		arr_g=vtk.vtkFloatArray()
		arr_g.SetName('G')

		arr_s=vtk.vtkFloatArray()
		arr_s.SetName('S')

		ipt=0
		for i in pg.panels:
			e=pg.panels[i]
			if not e.isWake or plotWake:
				pts=uniqueRows(e.pts)
				npts=len(pts)

				if npts ==4:
					shell=vtk.vtkQuad()
					pts=e.pts
				if npts==3:
					shell=vtk.vtkTriangle();
					#print e.gid
				if npts==3 or npts==4:
					for j in range(0,npts):
						ids=points.InsertNextPoint(pts[j])
					
						shell.GetPointIds().SetId(j,ipt)
						ipt+=1
					cells.InsertNextCell(shell)
					arr_eids.InsertNextValue(e.gid)
					arr_cp.InsertNextValue(e.pressureCoef)
					arr_g.InsertNextValue(e.G)
					arr_s.InsertNextValue(e.S)
			
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(points)
		polydata.SetPolys(cells)
		#polydata.Modified()

		
		polydata.GetCellData().AddArray(arr_eids)
		polydata.GetCellData().AddArray(arr_cp)
		#polydata.GetCellData().AddArray(arr_pload)
		polydata.GetCellData().AddArray(arr_s)
		polydata.GetCellData().AddArray(arr_g)
		return polydata


	mb = vtk.vtkMultiBlockDataSet() 
	mb.SetNumberOfBlocks(len(pg))
	kk=pg.keys()
	for i in range(0,len(kk)):
		mb.SetBlock( i, pgToPolyData(pg[kk[i]]) ) 
		mb.GetMetaData( i ).Set( vtk.vtkCompositeDataSet.NAME(), kk[i] ) 



	writer = vtk.vtkXMLMultiBlockDataWriter()
	writer.SetFileName(fname);
	writer.SetInputData(mb)
	writer.Write()
	del(writer)
	writer=None