import numpy
import math
class Protein:
    name=""
    baseline_deviation=[]
    baseline_mean=[]
    drug_values=[]
    cell_lines=[]
    drugs_applied=[]
    patient_drug_values=[]
    patient_baseline_mean=[]
    patient_baseline_deviation=[]
    patient_drugs_applied=[]
    patients=[]



drug_file="/Volumes/Samsung_T5 1/UDAI_PCA_BISSAN/LATEST_ANALYSIS/Files_From_Adam/SOC_postGAPDH_RawData_CellLines.txt"
baseline_file="/Volumes/Samsung_T5 1/UDAI_PCA_BISSAN/LATEST_ANALYSIS/Files_From_Adam/SOC_postGAPDH_RawBaseline_CellLines_FINAL.txt"

patient_baseline_file="/Volumes/Samsung_T5 1/UDAI_PCA_BISSAN/LATEST_ANALYSIS/Files_From_Adam/SOC_PatientData_RawBaseline.txt"
patient_drug_file="/Volumes/Samsung_T5 1/UDAI_PCA_BISSAN/LATEST_ANALYSIS/Files_From_Adam/SOC_PatientData_GAPDHnormalised.txt"

file=open(baseline_file,"r")

baseline_lines=file.readlines()

proteinlist=[]

linecount=0
celllinelist=[]
patientlist=[]

first_line_columns=baseline_lines[0].split("\t")

for column in first_line_columns:
    if column!="Cellline" and column!="Drug":
        newprotein = Protein()
        newprotein.name=column.rstrip("\n")
        newprotein.baseline_mean=[]
        newprotein.baseline_deviation=[]
        newprotein.drug_values=[]
        newprotein.cell_lines=[]
        newprotein.drugs_applied=[]
        newprotein.patient_baseline_deviation=[]
        newprotein.patient_baseline_mean=[]
        newprotein.patient_drug_values=[]
        newprotein.patient_drugs_applied=[]
        newprotein.patients=[]

        proteinlist.append(newprotein)


for a in  range(1,len(baseline_lines),6):
    columns=baseline_lines[a].rstrip("\n").split("\t")
    columncounter=0
    for b in range(2,len(columns)):
        proteinlist[columncounter].cell_lines.append(columns[0])
        proteinlist[columncounter].baseline_mean.append(numpy.mean([float(baseline_lines[a].split("\t")[b]),float(baseline_lines[a+1].split("\t")[b]),float(baseline_lines[a+2].split("\t")[b])]))
        proteinlist[columncounter].baseline_deviation.append(numpy.std([float(baseline_lines[a].split("\t")[b]),float(baseline_lines[a+1].split("\t")[b]),float(baseline_lines[a+2].split("\t")[b])]))


        columncounter=columncounter+1

file=open(drug_file,"r")

drug_lines=file.readlines()

for a in range(1,len(drug_lines)):
    columncounter = 0
    for b in range(2,len(drug_lines[a].rstrip("\n").split("\t"))):
        proteinlist[columncounter].drug_values.append(drug_lines[a].rstrip("\n").split("\t")[b])
        proteinlist[columncounter].drugs_applied.append(drug_lines[a].rstrip("\n").split("\t")[1])

        columncounter = columncounter + 1


drugslist=[]

for cellline in proteinlist[0].cell_lines:
    celllinelist.append(cellline)

for drugs in proteinlist[0].drugs_applied:
    if drugs not in drugslist:
        drugslist.append(drugs)



file=open(patient_baseline_file,"r")
baseline_lines=file.readlines()

for a in  range(1,len(baseline_lines),3):
    columns=baseline_lines[a].rstrip("\n").split("\t")
    columncounter=0
    for b in range(2,len(columns)):
        proteinlist[columncounter].patients.append(columns[0])
        proteinlist[columncounter].patient_baseline_mean.append(numpy.mean([float(baseline_lines[a].split("\t")[b]),float(baseline_lines[a+1].split("\t")[b]),float(baseline_lines[a+2].split("\t")[b])]))
        proteinlist[columncounter].patient_baseline_deviation.append(numpy.std([float(baseline_lines[a].split("\t")[b]),float(baseline_lines[a+1].split("\t")[b]),float(baseline_lines[a+2].split("\t")[b])]))


        columncounter=columncounter+1





file=open(patient_drug_file,"r")
drug_lines=file.readlines()


for a in range(1,len(drug_lines)):
    columncounter = 0
    for b in range(2,len(drug_lines[a].split("\t"))):
        proteinlist[columncounter].patient_drug_values.append(drug_lines[a].rstrip("\n").split("\t")[b])
        proteinlist[columncounter].patient_drugs_applied.append(drug_lines[a].rstrip("\n").split("\t")[1])

        columncounter = columncounter + 1

patient_drugslist=[]
for drugs in proteinlist[0].patient_drugs_applied:
    if drugs not in patient_drugslist:
        patient_drugslist.append(drugs)


for patient in proteinlist[0].patients:
    patientlist.append(patient)


output_file=open("CellLine_Fold_changes_output.txt","w")
output_file2=open("CellLine_Significant_hits_output.txt","w")

output_file.write("CellLine"+"\t"+"Drug"+"\t")
output_file2.write("CellLine"+"\t"+"Drug"+"\t")


for protein in proteinlist:
    output_file.write(protein.name+"\t")
    output_file2.write(protein.name + "\t")


output_file.write("\n")
output_file2.write("\n")
cellinecounter=0
for cellline in celllinelist:
    for a in range(0,len(drugslist)):
        output_file.write(cellline + "\t")
        output_file2.write(cellline + "\t")
        output_file.write(drugslist[a])
        output_file2.write(drugslist[a])

        output_file.write("\t")
        output_file2.write("\t")

        for protein in proteinlist:
            if protein.cell_lines[cellinecounter]==cellline:
                output_file.write(str(math.log2(float(protein.drug_values[a + len(drugslist) * (cellinecounter)]) / float(protein.baseline_mean[cellinecounter]))))

                if float(protein.drug_values[a+len(drugslist)*(cellinecounter)])>(float(protein.baseline_mean[cellinecounter])+2*float(protein.baseline_deviation[cellinecounter])) or float(protein.drug_values[a+len(drugslist)*(cellinecounter)])<(float(protein.baseline_mean[cellinecounter])-2*float(protein.baseline_deviation[cellinecounter])) :
                    output_file2.write(str(math.log2(float(protein.drug_values[a + len(drugslist) * (cellinecounter)]) / float(protein.baseline_mean[cellinecounter]))))
                else:
                    output_file2.write("0")


                output_file.write("\t")
                output_file2.write("\t")

        output_file.write("\n")
        output_file2.write("\n")


    cellinecounter=cellinecounter+1


output_file.close()
output_file2.close()

output_file=open("Patient_Fold_changes_output.txt","w")
output_file2=open("Patient_Significant_hits_output.txt","w")

output_file.write("Patient"+"\t"+"Drug"+"\t")
output_file2.write("Patient"+"\t"+"Drug"+"\t")


for protein in proteinlist:
    output_file.write(protein.name+"\t")
    output_file2.write(protein.name + "\t")


output_file.write("\n")
output_file2.write("\n")
cellinecounter=0
for patient in patientlist:
    for a in range(0,len(patient_drugslist)):
        output_file.write(patient + "\t")
        output_file2.write(patient + "\t")
        output_file.write(patient_drugslist[a])
        output_file2.write(patient_drugslist[a])

        output_file.write("\t")
        output_file2.write("\t")

        for protein in proteinlist:
            #print(protein.name)
            if protein.patients[cellinecounter]==patient:
                output_file.write(str(math.log2(float(protein.patient_drug_values[a + len(patient_drugslist) * (cellinecounter)]) / float(protein.patient_baseline_mean[cellinecounter]))))

                if float(protein.patient_drug_values[a+len(patient_drugslist)*(cellinecounter)])>(float(protein.patient_baseline_mean[cellinecounter])+2*float(protein.patient_baseline_deviation[cellinecounter])) or float(protein.patient_drug_values[a+len(patient_drugslist)*(cellinecounter)])<(float(protein.patient_baseline_mean[cellinecounter])-2*float(protein.patient_baseline_deviation[cellinecounter])) :
                    output_file2.write(str(math.log2(float(protein.patient_drug_values[a + len(patient_drugslist) * (cellinecounter)]) / float(protein.patient_baseline_mean[cellinecounter]))))
                else:
                    output_file2.write("0")


                output_file.write("\t")
                output_file2.write("\t")

        output_file.write("\n")
        output_file2.write("\n")


    cellinecounter=cellinecounter+1


output_file.close()
output_file2.close()
