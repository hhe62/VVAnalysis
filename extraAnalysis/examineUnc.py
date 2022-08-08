import ROOT
import pdb
import numpy as np

folders = ["16ErrorLog_2022-08-05","17ErrorLog_2022-08-04","18ErrorLog_2022-08-04"]
varstr = "nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass3j Mass34j Mass4j MassFull Mass0jFull Mass1jFull Mass2jFull Mass3jFull Mass34jFull Mass4jFull"
vars = varstr.split(" ")

def combineYears(l16,l17,l18,w16,w17,w18):
    lcorr = [w16*x+w17*y+w18*z for (x,y,z) in zip(l16,l17,l18)]
    luncorr = [((w16*x)**2+(w17*y)**2+(w18*z)**2)**0.5 for (x,y,z) in zip(l16,l17,l18)]
    return lcorr,luncorr

def sumListAbs(l):
    al = [abs(x) for x in l]
    return sum(al)

def analyzeYear(var,foldername,froot=None):
    dict = {}
    area = 0.
    fname = foldername+"/ErrorInfo_%s.log"%var
    
    hvar = froot.Get("tot_%s_unf"%var)
    area1 = hvar.Integral(1,hvar.GetNbinsX())
    #print("area from hist:%s"%area1)

    fin = open(fname)
    for line in fin:
        if "Area" in line:
            area = float(line.strip().split("Area: ")[1])
            #print("area from text:%s"%area)
        if "Source Up" in line:
            ln= line.strip().replace("Source Up ","")
            sys = ln.split(":")[0]
            contstr = (ln.split(":")[1][1:-1]).split(",") 
            cont = [float(x) for x in contstr]
            dict[sys]={} #Up occurs before Dn, so initialize here
            dict[sys]["Up"] = cont

        if "Source Dn" in line:
            ln= line.strip().replace("Source Dn ","")
            sys = ln.split(":")[0]
            contstr = (ln.split(":")[1][1:-1]).split(",") 
            cont = [float(x) for x in contstr]
            dict[sys]["Dn"] = cont

        if "Source Stat unc" in line:
            ln= line.strip()
            sys = "stat"
            contstr = (ln.split(":")[1][1:-1]).split(",") 
            cont = [float(x) for x in contstr]
            dict["stat"] = cont
        
        if "Source pdf unc" in line:
            ln= line.strip()
            sys = "pdf"
            contstr = (ln.split(":")[1][1:-1]).split(",") 
            cont = [float(x) for x in contstr]
            dict["pdf"] = cont
    
    if area == 0.:
        area = area1
    return area,dict

totDic = {}
areas = {}
for var in vars:
    totDic[var] = {}
    areas[var] = []
    for fd in folders:
        year = fd[0:2]
        froot = ROOT.TFile("%s/%s.root"%(fd,year))
        areay,dicty = analyzeYear(var,fd,froot)
        areas[var].append(areay)
        totDic[var][year] = dicty
        froot.Close()

dicComb = {}
jes_list = []
years = ["16","17","18"]
for var in vars:
    
    totarea = sum(areas[var])
    w16 = areas[var][0]/totarea
    w17 = areas[var][1]/totarea
    w18 = areas[var][2]/totarea
    #pdb.set_trace()
    fn_sys = []
    fn_corr = []
    fn_uncorr = []
    var_jes = 0.
    for sys in totDic[var]["18"].keys():

        if sys == "stat" or sys == "pdf":
            up16 = totDic[var]["16"][sys]
            up17 = totDic[var]["17"][sys]
            up18 = totDic[var]["18"][sys]
            upcorr,upuncorr = combineYears(up16,up17,up18,w16,w17,w18)
            unc_corr = sumListAbs(upcorr)
            unc_uncorr = sumListAbs(upuncorr)

        else:
            up16 = totDic[var]["16"][sys]["Up"]
            up17 = totDic[var]["17"][sys]["Up"]
            up18 = totDic[var]["18"][sys]["Up"]
            upcorr,upuncorr = combineYears(up16,up17,up18,w16,w17,w18)

            dn16 = totDic[var]["16"][sys]["Dn"]
            dn17 = totDic[var]["17"][sys]["Dn"]
            dn18 = totDic[var]["18"][sys]["Dn"]
            dncorr,dnuncorr = combineYears(dn16,dn17,dn18,w16,w17,w18)

            unc_corr = max(sumListAbs(upcorr),sumListAbs(dncorr))
            unc_uncorr = max(sumListAbs(upuncorr),sumListAbs(dnuncorr))

        fn_sys.append(sys)
        fn_corr.append(unc_corr)
        fn_uncorr.append(unc_uncorr)
        if sys == "jes":
            var_jes = unc_corr

    
    fn_corra = np.array(fn_corr)
    ind = (-fn_corra).argsort()
    fn_sys,fn_corr,fn_uncorr = [np.take(x,ind) for x in [fn_sys,fn_corr,fn_uncorr]]
    dicComb[var] = [fn_sys,fn_corr,fn_uncorr]
    jes_list.append(var_jes)

jes_lista = np.array(jes_list)
indjes = (-jes_lista).argsort()
vars_sort = vars
vars_sort = np.take(vars_sort,indjes)

for var in vars_sort:

    print("====%s==="%var)
    fn_sys,fn_corr,fn_uncorr = dicComb[var]
    for i in range(0,len(fn_sys)):
        print("%s %.4f %.4f"%(fn_sys[i],fn_corr[i],fn_uncorr[i]))





        
    