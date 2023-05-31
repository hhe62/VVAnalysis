import os
import sys
import json

#Usage: python Check_skimLogAndRedo.py folder year json_file

common_err = "WARNING: In non-interactive mode release checks e.g. deprecated releases, production architectures are disabled."
#Dir must be the root /nfs_scratch/<user>/<farmout_job_folder> for correct identification of folders
dir = sys.argv[1]
redolist = []
for root,dirs,files in os.walk(dir):
    for name in files:
        if "skim-" in name and ".err" in name:
            path = os.path.join(root,name)
            text = []
            with open(path) as ferr:
                for line in ferr:
                    text.append(line.rstrip())
            if (not text[0] in common_err) or len(text) >1:
                print("\nsomething happening with %s"%path)
                redolist.append(path.split("/")[1])
                print("Error Log:")
                print("")
                for item in text:
                    print(item)

#Redo the whole dataset if one file has issue
redolist = list(set(redolist))
print("Redo needed:")
print(redolist)

year = sys.argv[2]
#should end with .json
json_name = sys.argv[3]
rm_name = "cleanFailed_"+json_name.replace(".json",".sh")

with open("/hdfs/store/user/hhe62/%s"%rm_name,"w") as frm:
    for item in redolist:
        frm.write("rm -r %s\n"%item)

os.system("chmod u+x /hdfs/store/user/hhe62/%s"%rm_name)
print("/hdfs/store/user/hhe62/%s created"%rm_name)

with open("/hdfs/store/user/hehe/%s"%json_name) as json_file:
    obj = json.load(json_file)

#remove dataset that doesn't need resubmit from josn file 
match = False
for key in obj.keys():
    match = False
    for item in redolist:
        if key in item:
            match = True
    
    if not match:
        del obj[key]

with open("/hdfs/store/user/hehe/Resubmit_%s"%json_name,'w') as output_file:
    json.dump(obj,output_file,indent=4)

os.system("mv /hdfs/store/user/hehe/Resubmit_%s ~/vvanalysis_skim/CMSSW_10_3_1/src/Data_manager/ZZ4lRun2DatasetManager/FileInfo/ZZ4l%s/ntuples.json"%(json_name,year))
print("New json moved to ~/vvanalysis_skim/CMSSW_10_3_1/src/Data_manager/ZZ4lRun2DatasetManager/FileInfo/ZZ4l%s/ntuples.json"%year)
        



