import os, sys, json

dir = sys.argv[1]
redolist = []
for root,dirs,files in os.walk(dir):
    for name in files:
        if "skim-" in name and ".err" in name:
            path = os.path.join(root,name)
            text = []
            with open(path) as ferr:
                for line in ferr:
                    text.append(line.strip())
            if len(text) >1:
                print("\nsomething happening with %s"%path)
                redolist.append(path.split("/")[1])
                for item in text:
                    print(item)

redolist = list(set(redolist))
print("Redo needed:")
print(redolist)

year = sys.argv[2]
with open("/hdfs/store/user/hhe62/rmfailedfolders.sh","w") as frm:
    for item in redolist:
        frm.write("rm -r %s/\n"%item)

os.system("chmod u+x /hdfs/store/user/hhe62/rmfailedfolders.sh")
print("/hdfs/store/user/hhe62/rmfailedfolders.sh created")

with open("/hdfs/store/user/hehe/newNtuples%s.json"%year) as json_file:
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

with open("/hdfs/store/user/hehe/newNtuples%sResubmit.json"%year,'w') as output_file:
        json.dump(obj,output_file,indent=4)

os.system("mv /hdfs/store/user/hehe/newNtuples%sResubmit.json ~/vvanalysis_skim/CMSSW_10_3_1/src/Data_manager/ZZ4lRun2DatasetManager/FileInfo/ZZ4l%s/ntuples.json"%(year,year))
print("New json moved to ~/vvanalysis_skim/CMSSW_10_3_1/src/Data_manager/ZZ4lRun2DatasetManager/FileInfo/ZZ4l%s/ntuples.json"%year)
        



