import os
import sys

#Usage: python Check_skimLogAndRedo.py folder 

common_err = "WARNING: In non-interactive mode release checks e.g. deprecated releases, production architectures are disabled."
#dir must be the root /nfs_scratch/<user>/<farmout_job_folder> for correct identification of folders
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
