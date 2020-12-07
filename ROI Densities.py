#modifieid the roi_density.py script assumes hu to density table uses linear interpolation

# Script to compute the minimum, maximum and average density and the standard deviation
# of the density for an ROI for the currently open patient
# Select the ROI to compute average density for by setting the parameter 'roi_name'
# Density override values are not taken into account in the computation

from connect import *
import sys
from time import time

connect(22392)

#roi_name = 'CTV'
def FindImagingSystem(frameofreference):
    patient=get_current("Patient")
    for exam in patient.Examinations:
        if frameofreference==exam.EquipmentInfo.FrameOfReference:
            return exam.EquipmentInfo.ImagingSystemReference.ImagingSystemName,exam.Name
    return None
def interpolate(yn,x1,y1,x2,y2):
    #pt to interpolate to followed by coordinates
    if y2-y1==0:
        return 0
    else:
        return (x1+(yn-y1)*(x2-x1)/(y2-y1))
def getVolume(roiname,examname):
    patient=get_current("Patient")
    return patient.PatientModel.StructureSets[examname].RoiGeometries[roiname].GetRoiVolume()
def getCT(ctname):

    machine_db=get_current("MachineDB")
    for ct in machine_db.ImagingSystems:
        if ct.Name==ctname:
            return ct

def getHUtoRHO(ctname):
    ct=getCT(ctname)
    hu=ct.CtToDensityTable.HU
    rho=ct.CtToDensityTable.Rho
    return [hu,rho]





def CTtoHU(data,density):


    #data=getHUtoRHO(ctname)
    i=0

    while i <len(data[0])-1:
        #print data[0][i],data[1][i]
        if density >= data[1][i] and density<=data[1][i+1]:
            #print density,'converted',interpolate(density,data[0][i],data[1][i],data[0][i+1],data[1][i+1])
            return interpolate(density,data[0][i],data[1][i],data[0][i+1],data[1][i+1])



        i+=1




def GetHUData(roi_name):
    try:
        beam_set = get_current('BeamSet')
    except:
        print 'No beam set is currently loaded. Cannot compute average density.'
        sys.exit()

    dgr = beam_set.FractionDose.GetDoseGridRoi(RoiName=roi_name)
    if not dgr:
        print 'No ROI named {0} exists.'.format(roi_name)
        sys.exit()

    if not beam_set.FractionDose.OnDensity.DensityValues:
        print 'No density values are computed for the current dose grid. Cannot compute average density.'
        sys.exit()

    ctname,examname=FindImagingSystem(beam_set.FrameOfReference)
    if ctname==None:
        print 'No reference CT found.'
        sys.exit()
    roivol=getVolume(roi_name,examname)
    ctdata=getHUtoRHO(ctname)
    density = [d for d in beam_set.FractionDose.OnDensity.DensityValues.DensityData]
    average_density = 0
    min_density = 25
    max_density = 0
    avg_hu=0
    min_hu=3000.0
    max_hu=-3000.0
    #print 'start here'
    quickvar=[]
    for i, v in zip(dgr.RoiVolumeDistribution.VoxelIndices, dgr.RoiVolumeDistribution.RelativeVolumes):
        den=density[i]
        #average_density += den * v#sity[i] * v
        hu=CTtoHU(ctdata,den)#sity[i])
        avg_hu+=hu*v#CTtoHU(ctname,density[i]) * v
        #min_density = min(min_density, den)#density[i])
        #max_density = max(max_density, den)#density[i])
        #if roi_name=='bladder':
        #    print min_hu,hu
        min_hu = min(min_hu,hu)# CTtoHU(ctname,density[i]))
        max_hu = max(max_hu,hu)# CTtoHU(ctname,density[i]))
        #print average_density,avg_hu
        quickvar.append([hu,v])

    #print 'Looking for variance here'
    huvariance = 0
    variance = 0
    for line in quickvar:
        huvariance+=line[1]*(line[0]-avg_hu)**2
#    for i, v in zip(dgr.RoiVolumeDistribution.VoxelIndices, dgr.RoiVolumeDistribution.RelativeVolumes):
#        #variance += v * (density[i] - average_density)**2
#        huvariance+=v*(CTtoHU(ctname,density[i]) - avg_hu)**2


    standard_deviation = variance ** 0.5
    stdevhu=huvariance ** 0.5

    #print 'Mean density of ROI {0}: {1} g / cm^3'.format(roi_name, average_density)
    #print 'Standard deviation of density for ROI {0}: {1} g / cm^3'.format(roi_name, standard_deviation)
    #print 'Minimum density of ROI {0}: {1} g / cm^3'.format(roi_name, min_density)
    #print 'Maximum density of ROI {0}: {1} g / cm^3'.format(roi_name, max_density)
    #return average_density, standard_deviation, min_density, max_density,avg_hu,stdevhu,min_hu,max_hu
    return roivol,avg_hu,stdevhu,min_hu,max_hu

#this one did not work out.. binary data seems to have every other one a 0, and not sure that it's actual HU
def oldGetHUData(roi_name):
    try:
        beam_set = get_current('BeamSet')
    except:
        print 'No beam set is currently loaded. Cannot compute average density.'
        sys.exit()

    dgr = beam_set.FractionDose.GetDoseGridRoi(RoiName=roi_name)
    if not dgr:
        print 'No ROI named {0} exists.'.format(roi_name)
        sys.exit()

    if not beam_set.FractionDose.OnDensity.DensityValues:
        print 'No density values are computed for the current dose grid. Cannot compute average density.'
        sys.exit()

    patient=get_current("Patient")
    #density = [d for d in beam_set.FractionDose.OnDensity.DensityValues.DensityData]
    density = [d for d in patient.Examinations[0].Series[0].ImageStack.PixelData]
    average_density = 0
    min_density = 4000
    max_density = -2000
    data=[]
    for i, v in zip(dgr.RoiVolumeDistribution.VoxelIndices, dgr.RoiVolumeDistribution.RelativeVolumes):
        average_density += density[i] * v
        #print density[i],v
#        min_density = min(min_density, density[i])
#        max_density = max(max_density, density[i])
        data.append(density[i])

    variance = 0
    for i, v in zip(dgr.RoiVolumeDistribution.VoxelIndices, dgr.RoiVolumeDistribution.RelativeVolumes):
        variance += v * (density[i] - average_density)**2

    standard_deviation = variance ** 0.5
    print

    print 'Mean density of ROI {0}: {1} g / cm^3'.format(roi_name, average_density)
    print 'Standard deviation of density for ROI {0}: {1} g / cm^3'.format(roi_name, standard_deviation)
    print 'Minimum density of ROI {0}: {1} g / cm^3'.format(roi_name, min_density)
    print 'Maximum density of ROI {0}: {1} g / cm^3'.format(roi_name, max_density)
    return average_density, standard_deviation, min_density, max_density

t0=time()
patient = get_current("Patient")
#beam_set=get_current("BeamSet")

roistats=['Organ Name','Organ Type','Volume','Average HU','Standard Dev HU','Minimum HU','Maximum HU','Analysis Time']
#name='PROSTATE'

#avg,stdev,min,max,avg_hu,stdevhu,min_hu,max_hu=GetDensityData(name)
#roistats.append([name,avg,stdev,min,max,avg_hu,stdevhu,min_hu,max_hu])
#roistats.append(['Name','density average',])
#print roistats[0]
#crashnow
last=time()
for roi in patient.PatientModel.RegionsOfInterest:
    try:
        #avg,stdev,min,max,
        roivol,avg_hu,stdevhu,min_hu,max_hu=GetHUData(roi.Name)
        roistats.append([roi.Name,roi.OrganData.OrganType,roivol,avg_hu,stdevhu,min_hu,max_hu,time()-last])
        #print roistats[-1]#[roi.Name,avg_hu,stdevhu,min_hu,max_hu]
    except:
        roistats.append([roi.Name,'Failed'])
    last=time()
    #print "\nProcessing Time: "+str(last)+" seconds"
for line in roistats:
    print line
print "\nProcessing Time: "+str(time()-t0)+" seconds"