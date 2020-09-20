import vtk
import math
import os
import numpy as np

def test_file_writer_output(filename):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)

    reader.Update()

    pdata = reader.GetOutput()

    assert pdata.GetNumberOfPoints() == 11137

    for i in range(pdata.GetNumberOfPoints()):
        if(pdata.GetPointData().GetArray('Boundary').GetValue(i) == 0):

            vx = pdata.GetPointData().GetArray('Velocity').GetTuple(i)[0]
            vy = pdata.GetPointData().GetArray('Velocity').GetTuple(i)[1]
            assert math.sqrt(vx * vx + vy * vy) < 20 

            assert (pdata.GetPointData().GetArray('Pressure').GetValue(i)< 200000)
            p_x = pdata.GetPoint(i)[0]
            p_y = pdata.GetPoint(i)[1]
            assert ((p_x >= 0) and (p_x <= 20))
            # assert ((p_y >= 0) and (p_y <= 10))


for j in range(11137):
        path = "./data_"+ str(j) +".vtp" 
        test_file_writer_output(path)
        print("step"+ str(j)+"Pass!")

