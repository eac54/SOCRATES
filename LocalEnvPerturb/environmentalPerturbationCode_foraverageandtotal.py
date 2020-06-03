import os

class Protein:
    def __init__(self):
        self.name=""
        self.interactors=[]
        self.phosphochange=0
        self.envpertscore=0



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

celllinenamelist=[]
celllines=[]

for file in os.listdir(input_directory_name):

    if file.endswith("txt")  and file.startswith(".")==False:
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
        
        newprotein=Protein()
        newprotein.name="BRAF"
        newprotein.phosphochange=0

        temp_proteinlist.append(newprotein)
        temp_proteinlistname.append(newprotein.name)

        newprotein=Protein()
        newprotein.name="HSP90"
        newprotein.phosphochange=0

        temp_proteinlist.append(newprotein)
        temp_proteinlistname.append(newprotein.name)

        newprotein=Protein()
        newprotein.name="PIK3CA"
        newprotein.phosphochange=0

        temp_proteinlist.append(newprotein)
        temp_proteinlistname.append(newprotein.name)

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
            protein.envpertscore=abs(protein.phosphochange)
            for interactor in protein.interactors:
                protein.envpertscore=protein.envpertscore+abs(interactor.phosphochange)

            protein.envpertscore=protein.envpertscore/len(protein.interactors) #ACTIVATE FOR AVERAGE CALCULATION

outfilelocation=output_directory_name+"perturbation_matrix_averageperturbation.txt"
output=open(outfilelocation,"w")

output.write("\t")

for protein in celllines[0].proteinlist[0]:
    output.write(protein.name+"\t")

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



