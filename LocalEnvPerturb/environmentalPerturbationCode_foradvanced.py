import os
import numpy as np

class Protein:
    def __init__(self):
        self.name=""
        self.interactors=[]
        self.phosphochange=0
        self.envpertscore=0
        self.degree=0



class Cellline:
    def __init__(self):
        self.name=""
        self.drugapplied=[]
        self.proteinlist=[]



input_directory_name="/Volumes/Samsung_T5/SOCRATES/Advanced_analysis_files/EnvironmentalPerturbationCode/input_final/"
output_directory_name="/Volumes/Samsung_T5/SOCRATES/Advanced_analysis_files/EnvironmentalPerturbationCode/output_final/"

interactome_file="/Volumes/Samsung_T5/SOCRATES/Advanced_analysis_files/EnvironmentalPerturbationCode/2018Interactome_FINAL_UdaiNetwork_LizziesInteractionsAdded.txt"

file=open(interactome_file,"r")

interactome_lines=file.readlines()
proteinlist=[]
proteinlistname=[]

for line in interactome_lines:
    columns=line.rstrip("\n").split("\t")
    if columns[0]!=columns[2]:
        if columns[0] not in proteinlistname:
            newprotein=Protein()

            newprotein.name=columns[0]
            newprotein.interactors.append(columns[2])

            proteinlist.append(newprotein)
            proteinlistname.append(newprotein.name)
        else:
            protindex=proteinlistname.index(columns[0])
            if columns[2] not in proteinlist[protindex].interactors:
                proteinlist[protindex].interactors.append(columns[2])

        if columns[2] not in proteinlistname:
            newprotein = Protein()

            newprotein.name = columns[2]
            newprotein.interactors.append(columns[0])

            proteinlist.append(newprotein)
            proteinlistname.append(newprotein.name)
        else:
            protindex = proteinlistname.index(columns[2])
            if columns[0] not in proteinlist[protindex].interactors:
                proteinlist[protindex].interactors.append(columns[0])

distance_matrix=np.empty((len(proteinlistname),len(proteinlistname)))

for protein in proteinlist:
    tempnumberlist=[]

    rowindex=proteinlistname.index(protein.name)

    distance_matrix[rowindex, rowindex]=-1
    tempnumberlist.append(rowindex)

    for interactor in protein.interactors:
        columnindex=proteinlistname.index(interactor)
        distance_matrix[rowindex,columnindex]=0
        tempnumberlist.append(columnindex)


    for interactor in protein.interactors:
        columnindex = proteinlistname.index(interactor)
        for scnddegree in proteinlist[columnindex].interactors:
            scnddegreeindex=proteinlistname.index(scnddegree)
            if proteinlistname[rowindex]!=proteinlistname[scnddegreeindex] and scnddegreeindex not in tempnumberlist:
                distance_matrix[rowindex, scnddegreeindex] = 1
            if scnddegreeindex not in tempnumberlist:
                tempnumberlist.append(scnddegreeindex)

    for interactor in protein.interactors:
        columnindex = proteinlistname.index(interactor)
        for scnddegree in proteinlist[columnindex].interactors:
            scnddegreeindex=proteinlistname.index(scnddegree)
            for thirddegree in proteinlist[scnddegreeindex].interactors:
                thirddegreeindex = proteinlistname.index(thirddegree)

                if proteinlistname[rowindex] != proteinlistname[thirddegreeindex] and thirddegreeindex not in tempnumberlist:
                    distance_matrix[rowindex, thirddegreeindex] = 2
                if thirddegreeindex not in tempnumberlist:
                    tempnumberlist.append(thirddegreeindex)



celllinenamelist=[]
celllines=[]

for file in os.listdir(input_directory_name):

    if file.endswith("txt") and file.startswith(".")==False:
        filename=input_directory_name+file
        phosphodata=open(filename,"r")
        Phospho_lines=phosphodata.readlines()

        temp_proteinlist=[]
        temp_proteinlistname = []
        linenumber = 0

        for line in Phospho_lines:
            if linenumber!=0:
                columns=line.rstrip("\n").split("\t")


                newprotein=Protein()
                newprotein.name=columns[0]
                newprotein.phosphochange=float(columns[1])

                temp_proteinlist.append(newprotein)
                temp_proteinlistname.append(newprotein.name)

            linenumber=linenumber+1
            print(linenumber)

        for protein in temp_proteinlist:
            protein_index = proteinlistname.index(protein.name)

            for interactor in proteinlist[protein_index].interactors:
                if interactor not in ["ERK2","GAPDH"]:#Because not in the phospholist
                    temp_index=temp_proteinlistname.index(interactor)
                    protein.interactors.append(temp_proteinlist[temp_index])



        if file.split("_")[0] not in celllinenamelist:
            newcellline=Cellline()

            newcellline.name=file.split("_")[0]
            newcellline.drugapplied.append(file.split("_")[1].split(".")[0])
            newcellline.proteinlist.append(temp_proteinlist)

            celllines.append(newcellline)
            celllinenamelist.append(newcellline.name)
        else:
            cellline_index=celllinenamelist.index(file.split("_")[0])

            celllines[cellline_index].drugapplied.append(file.split("_")[1].split(".")[0])

            celllines[cellline_index].proteinlist.append(temp_proteinlist)




letmeknow=True

for cellline in celllines:
    for a in range(0,len(cellline.drugapplied)):
        for protein in cellline.proteinlist[a]:
            firstindex=proteinlistname.index(protein.name)
            protein.envpertscore=protein.envpertscore+abs(protein.phosphochange)*np.exp(0)
            for scndprotein in cellline.proteinlist[a]:
                secondindex = proteinlistname.index(scndprotein.name)
                if -1<distance_matrix[firstindex,secondindex]<5:
                    protein.envpertscore = protein.envpertscore + abs(scndprotein.phosphochange) * np.exp(-distance_matrix[firstindex,secondindex])

outfilelocation=output_directory_name+"perturbation_matrix_MOSTDETAILEDPERTURBATION.txt"
output=open(outfilelocation,"w")
output2=open(output_directory_name+"_proteindegrees.txt","w")

output.write("\t")

for protein in celllines[0].proteinlist[0]:
    output.write(protein.name+"\t")
    output2.write(protein.name+"\t")
    output2.write(str(len(protein.interactors))+"\n")

output.write("\n")


for cellline in celllines:
    for a in range(0,len(cellline.drugapplied)):
        print(str(a))
        output.write(cellline.name +"_"+cellline.drugapplied[a]+ "\t")
        for protein in cellline.proteinlist[a]:
            output.write(str(protein.envpertscore))
            output.write("\t")

        output.write("\n")



output.close()



