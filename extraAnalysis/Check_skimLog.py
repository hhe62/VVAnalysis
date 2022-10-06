import os,sys

dir = sys.argv[1]
for root,dirs,files in os.walk(dir):
    for name in files:
        if "skim-" in name and ".err" in name:
            path = os.path.join(root,name)
            text = []
            with open(path) as ferr:
                for line in ferr:
                    text.append(line.strip())
            if len(text) >1:
                print("something happening with %s"%path)
                for item in text:
                    print(item)