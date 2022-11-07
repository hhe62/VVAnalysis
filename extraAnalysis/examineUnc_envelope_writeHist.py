import ROOT
import pdb,json
import numpy as np
import array

#======================================
#Name folders starting with 16,17,18, put 16,17,18.root into folders
#======================================
folders = ["16Error","17Error","18Error"]
varstr = "nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass3j Mass34j Mass4j MassFull Mass0jFull Mass1jFull Mass2jFull Mass3jFull Mass34jFull Mass4jFull"
vars = varstr.split(" ")
vars = ["nJets"]

def combineYears(l16,l17,l18,w16,w17,w18):
    lcorr = [w16*x+w17*y+w18*z for (x,y,z) in zip(l16,l17,l18)]
    luncorr = [((w16*x)**2+(w17*y)**2+(w18*z)**2)**0.5 for (x,y,z) in zip(l16,l17,l18)]
    return lcorr,luncorr

def sumListAbs(l):
    al = [abs(x) for x in l]
    return sum(al)

def sqrt_sum(l1,l2):
    l = [(x**2 + y**2)**0.5 for (x,y) in zip(l1,l2)]
    return l

def sep_up_dn(lu,ld):
    lun = []
    ldn = []
    for x,y in zip(lu,ld):
        lun.append(max(x,y))
        ldn.append(min(x,y))
    
    return lun,ldn

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
FillDic = {} #For storing tot up and down unc bin content
areas = {}
for var in vars:
    totDic[var] = {}
    areas[var] = []
    FillDic[var] = []
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
    tot_corrUp,tot_uncorrUp = [],[]
    tot_corrDn,tot_uncorrDn = [],[]
    for sys in totDic[var]["18"].keys():

        if sys == "stat" or sys == "pdf":
            up16 = totDic[var]["16"][sys]
            up17 = totDic[var]["17"][sys]
            up18 = totDic[var]["18"][sys]
            upcorr,upuncorr = combineYears(up16,up17,up18,w16,w17,w18)
            unc_corr = sumListAbs(upcorr)
            unc_uncorr = sumListAbs(upuncorr)

            if tot_corrUp == []:
                if sys == "stat":
                    tot_corrUp,tot_uncorrUp = upuncorr,upuncorr
                    tot_corrDn,tot_uncorrDn = upuncorr,upuncorr
                else:
                    tot_corrUp,tot_uncorrUp = upcorr,upcorr
                    tot_corrDn,tot_uncorrDn = upcorr,upcorr

            else:
                if sys == "stat":
                    tot_corrUp,tot_uncorrUp = sqrt_sum(tot_corrUp,upuncorr),sqrt_sum(tot_uncorrUp,upuncorr)
                    tot_corrDn,tot_uncorrDn = sqrt_sum(tot_corrDn,upuncorr),sqrt_sum(tot_uncorrDn,upuncorr)
                else:
                    tot_corrUp,tot_uncorrUp = sqrt_sum(tot_corrUp,upcorr),sqrt_sum(tot_uncorrUp,upcorr)
                    tot_corrDn,tot_uncorrDn = sqrt_sum(tot_corrDn,upcorr),sqrt_sum(tot_uncorrDn,upcorr)



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

            #change up/down properly for tot unc. calculation
            upcorr,dncorr = sep_up_dn(upcorr,dncorr)
            upuncorr, dnuncorr = sep_up_dn(upuncorr,dnuncorr)

            if tot_corrUp == []:
                if sys == "jes":
                    tot_corrUp,tot_uncorrUp = upcorr,upuncorr
                    tot_corrDn,tot_uncorrDn = dncorr,dnuncorr
                if sys == "jer":
                    tot_corrUp,tot_uncorrUp = upuncorr,upuncorr
                    tot_corrDn,tot_uncorrDn = dnuncorr,dnuncorr
                else:
                    tot_corrUp,tot_uncorrUp = upcorr,upcorr
                    tot_corrDn,tot_uncorrDn = dncorr,dncorr

            else:
                if sys == "jes":
                    tot_corrUp,tot_uncorrUp = sqrt_sum(tot_corrUp,upcorr),sqrt_sum(tot_uncorrUp,upuncorr)
                    tot_corrDn,tot_uncorrDn = sqrt_sum(tot_corrDn,dncorr),sqrt_sum(tot_uncorrDn,dnuncorr)

                if sys == "jer":
                    tot_corrUp,tot_uncorrUp = sqrt_sum(tot_corrUp,upuncorr),sqrt_sum(tot_uncorrUp,upuncorr)
                    tot_corrDn,tot_uncorrDn = sqrt_sum(tot_corrDn,dnuncorr),sqrt_sum(tot_uncorrDn,dnuncorr)

                else:
                    tot_corrUp,tot_uncorrUp = sqrt_sum(tot_corrUp,upcorr),sqrt_sum(tot_uncorrUp,upcorr)
                    tot_corrDn,tot_uncorrDn = sqrt_sum(tot_corrDn,dncorr),sqrt_sum(tot_uncorrDn,dncorr)

        fn_sys.append(sys)
        fn_corr.append(unc_corr)
        fn_uncorr.append(unc_uncorr)
        if sys == "jes":
            var_jes = unc_corr

    totunc_corr = max(sumListAbs(tot_corrUp),sumListAbs(tot_corrDn))
    totunc_uncorr = max(sumListAbs(tot_uncorrUp),sumListAbs(tot_uncorrDn))

    fn_corra = np.array(fn_corr)
    ind = (-fn_corra).argsort()
    fn_sys,fn_corr,fn_uncorr = [np.take(x,ind) for x in [fn_sys,fn_corr,fn_uncorr]]
    dicComb[var] = [fn_sys,fn_corr,fn_uncorr,totunc_corr,totunc_uncorr]
    #jes_list is used to append the quantity used to sort variables
    jes_list.append(abs(totunc_corr-totunc_uncorr)/totunc_corr)
    #jes_list.append(var_jes)

    FillDic[var] = [tot_corrUp, tot_corrDn]

jes_lista = np.array(jes_list)
indjes = (-jes_lista).argsort()
vars_sort = vars
vars_sort = np.take(vars_sort,indjes)

for var in vars_sort:

    print("====%s==="%var)
    print("%-10s %-6s uncorr"%(" ","corr"))
    fn_sys,fn_corr,fn_uncorr,final_corr,final_uncorr = dicComb[var]
    for i in range(0,len(fn_sys)):
        print("%-10s %.4f %.4f"%(fn_sys[i],fn_corr[i],fn_uncorr[i]))
    
    print("Total uncertainty with jes correlated:%.4f uncorrelated:%.4f, relative diff %.4f"%(final_corr,final_uncorr, abs(final_corr-final_uncorr)/final_corr))

with open('varsFile.json') as var_json_file:
    myvar_dict = json.load(var_json_file)

#Files with tot unc hists removed
fr2 = ROOT.TFile("Full.root","READ")
fout = ROOT.TFile("out.root","UPDATE")
for var in vars:

    _binning = myvar_dict[var]["_binning"]
    histbins=array.array('d',_binning)
    hUncUp=ROOT.TH1D("tot_%s_totUncUp"%var,"Total Up Uncert.",len(histbins)-1,histbins)
    hUncDn=ROOT.TH1D("tot_%s_totUncDown"%var,"Total Dn Uncert.",len(histbins)-1,histbins)
    for i in range(1,hUncUp.GetNbinsX()+1):
        totUncUp = FillDic[var][0][i-1]
        totUncDn = FillDic[var][1][i-1]
        hUncUp.SetBinContent(i,totUncUp*totarea)
        hUncDn.SetBinContent(i,totUncDn*totarea)
    
    #Command line tool doesn't work for some hist so have to do it manually
    fr2.cd()
    olddata = fr2.Get("tot_%s_data"%var)
    oldtrue = fr2.Get("tot_%s_true"%var)
    oldtrueAlt = fr2.Get("tot_%s_trueAlt"%var)
    oldBkg = fr2.Get("tot_%s_bkg"%var)
    oldUnf = fr2.Get("tot_%s_unf"%var)
    olddSigMC= fr2.Get("tot_%s_SigMC"%var)

    fout.cd()
    tmpHists = [olddata,oldtrue,oldtrueAlt,oldBkg,oldUnf,olddSigMC,hUncUp,hUncDn]
    for h in tmpHists:
        h.Write()
  

fr2.Close()



        
    
