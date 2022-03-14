#Partially reimplement RooUnfold package methods in python
import ROOT
import math,pdb

class RooUnfoldResponse(object):

    def __init__(self,measured, truth, response):
        self.Init()
        self.SetupMatrix (measured, truth, response)

    def Init(self):
        self._overflow= 0
        self.Setup()

    def Setup(self):
        self._tru= self._mes= self._fak= 0 # not touched in clearCache()
        self._res= 0             # not touched in clearCache()
        self._vMes= self._eMes= self._vFak= self._vTru= self._eTru= 0
        self._mRes= self._eRes= 0
        self._nm= self._nt= self._mdim= self._tdim= 0 # not touched in clearCache()
        self._cached= False

    def ClearCache(self):
        self._vMes= 0
        self._eMes= 0
        self._vFak= 0
        self._vTru= 0
        self._eTru= 0
        self._mRes= 0
        self._eRes= 0
        self._cached= False

    def Reset(self):

    # Resets object to initial state.
        self.ClearCache()
        del self._mes
        del self._fak
        del self._tru
        del self._res
        self.Setup()

    @staticmethod
    def GetBin (h, i, overflow):

    # vector index (0..nx*ny-1) -> multi-dimensional histogram
    # global bin number (0..(nx+2)*(ny+2)-1) skipping under/overflow bins
        return ( i+(0 if overflow else 1) if h.GetDimension()<2 else Exception("hist dimension problem"))

    @staticmethod
    def H2V(h, nb, overflow):
    #Returns TVectorD of the bin contents of the input histogram
        if (overflow):
            nb += 2
        v= ROOT.TVectorD(nb)
        if (not h): 
            return v
        for i in range(0,nb):
            v[i]= h.GetBinContent( RooUnfoldResponse.GetBin(h, i, overflow))
  
        return v
    
    def Vtruth(self):
    #Truth distribution as a TVectorD
        if (not self._vTru):
            self._vTru= RooUnfoldResponse.H2V(self._tru, self._nt, self._overflow)
            self._cached= True if self._vTru else False
        return self._vTru
    
    def Vfakes(self):
    #Fakes distribution as a TVectorD
        if (not self._vFak):
            self._vFak= RooUnfoldResponse.H2V(self._fak, self._nm, self._overflow)
            self._cached= True if self._vFak else False 
        return self._vFak

    def GetNbinsMeasured(self):
    # Total number of bins in the measured distribution
        return self._nm
    
    def GetNbinsTruth(self):
    #Total number of bins in the truth distribution
        return self._nt

    def FakeEntries(self):
    #Return number of fake entries
        return self._fak.GetEntries() if self._fak else 0.0

    def UseOverflowStatus(self):
    #Get UseOverflow setting
        return self._overflow

    def Hfakes(self):
        return self._fak

    def Htruth(self):
    #Truth distribution, used for normalisation
        return self._tru
    
    def Hmeasured(self):
    #Measured distribution, including fakes ->? isn't fake stored in separate _fak?
        return self._mes

    def Hresponse(self):
    #Response matrix as a 2D-histogram: (x,y)=(measured,truth)
        return self._res



    def SetupMatrix (self, measured, truth, response):

    # Set up from already-filled histograms.
    # "response" gives the response matrix, measured X truth.
    # "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
    # but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
    # in "truth" for unmeasured events (inefficiency).
    # "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
    # to indicate, respectively, no fakes and/or no inefficiency.
        self.Reset()
        oldstat= ROOT.TH1.AddDirectoryStatus()
        ROOT.TH1.AddDirectory(False)
        self._res= response.Clone()
        if (measured): 
            self._mes=  measured.Clone()
            self._fak=  measured.Clone("fakes")
            self._fak.Reset()
            self._fak.SetTitle("Fakes")
            self._mdim= self._mes.GetDimension()
    
        if (truth): 
            self._tru=   truth.Clone()
            self._tdim= self._tru.GetDimension()
    
        ROOT.TH1.AddDirectory (oldstat)
        if (self._overflow and (self._mdim > 1 or self._tdim > 1)):
            print("UseOverflow setting ignored for multi-dimensional distributions") 
            self._overflow= 0
    
        self._nm= self._mes.GetNbinsX() * self._mes.GetNbinsY() * self._mes.GetNbinsZ()
        self._nt= self._tru.GetNbinsX() * self._tru.GetNbinsY() * self._tru.GetNbinsZ()
        if (self._nm != self._res.GetNbinsX() or self._nt != self._res.GetNbinsY()):
            print( "Warning: RooUnfoldResponse measured X truth is " ,self._nm ," X ", self._nt,
            ", but matrix is ",self._res.GetNbinsX()," X ",self._res.GetNbinsY())
            raise Exception('Something wrong in dimension')
    
        first=1
        nm= self._nm
        nt= self._nt
        s= self._res.GetSumw2N()

        if (self._overflow):
            first= 0
            nm += 2
            nt += 2
    

        if (not measured or self._mes.GetEntries() == 0.0):
            raise Exception("No measured hist content")
        else:
        #Fill fakes from the difference of self._mes - self._res.ProjectionX()
        #Always include under/overflows in sum of truth.
            sm= self._mes.GetSumw2N()
            nfake=0
            for i in range(0,nm):
                nmes= 0.0
                wmes= 0.0
                for j in range(0,self._nt+2):
                    nmes += self._res.GetBinContent (i+first, j)
                    if (s): 
                        wmes += pow (self._res.GetBinError(i+first, j), 2)
        
                bin= RooUnfoldResponse.GetBin(self._mes, i, self._overflow)
                fake= self._mes.GetBinContent (bin) - nmes
            #cout<<i+1<<":inside"<<self._mes.GetBinContent (bin)<<endl
                if (fake!=0.0): 
                    nfake+=1
                if (not s):
                    wmes= nmes
                self._fak.SetBinContent (bin, fake)
                self._fak.SetBinError(bin, math.sqrt(wmes + (pow(self._mes.GetBinError(bin),2) if sm else self._mes.GetBinContent(bin))))
        
            #if ROOT_VERSION_CODE >= ROOT_VERSION(5,13,0)
            self._fak.SetEntries(self._fak.GetEffectiveEntries())  # 0 entries if 0 fakes

    

        if (not truth or self._tru.GetEntries() == 0.0):
            raise Exception("No truth hist content")

class RooUnfold(object):
    def __init__(self, res, meas):
        self.Init()
        self.Setup (res, meas)

    def Init(self):
        self._res= self._resmine= 0
        self._vMes= self._eMes= 0
        self._covMes= self._covL= 0
        self._meas= self._measmine= 0
        self._nm= self._nt= 0
        self._verbose= 1
        self._overflow= 0
        self._dosys= self._unfolded= self._haveCov= self._haveCovMes= self._fail= self._have_err_mat= self._haveErrors= self._haveWgt= False
        self._withError= 'kDefault'
        self._NToys=50
        self.GetSettings()
    
    def GetSettings(self):
    #Gets maximum and minimum parameters and step size
        self._minparm=0
        self._maxparm=0
        self._stepsizeparm=0
        self._defaultparm=0

    def Setup (self,res, meas):
        self.Reset()
        self.SetResponse (res)
        self.SetMeasured (meas)

    def Reset(self):
        self.Destroy()
        self.Init()

    def Destroy(self):
        del self._measmine
        del self._vMes
        del self._eMes
        del self._covMes
        del self._covL
        del self._resmine
    
    def SetResponse (self,res):
    # Set response matrix for unfolding.
        self._resmine= 0
        self._res= res
        self._overflow= 1 if self._res.UseOverflowStatus() else 0
        self._nm= self._res.GetNbinsMeasured()
        self._nt= self._res.GetNbinsTruth()
        if (self._overflow):
            self._nm += 2
            self._nt += 2
    
    def SetMeasured (self,meas):
    #Set measured distribution and errors. RooUnfold does not own the histogram.
        self._meas= meas
        self._vMes= 0
        self._eMes= 0

    def Hreco(self,withError='kErrors'):
    #Creates reconstructed distribution. Error calculation varies by withError:
    #0: No errors
    #1: Errors from the square root of the diagonals of the covariance matrix given by the unfolding
    #2: Errors from the square root of of the covariance matrix given by the unfolding
    #3: Errors from the square root of the covariance matrix from the variation of the results in toy MC tests
    
        reco= self._res.Htruth().Clone('Bayes_reco')
        reco.Reset()
        reco.SetTitle('Bayes_reco_title')
        if (not self.UnfoldWithErrors (withError)):
            withError= 'kNoError'
        if (not self._unfolded):
            return reco

        for i in range(0,self._nt):
            j= RooUnfoldResponse.GetBin (reco, i, self._overflow)
            reco.SetBinContent (j, self._rec(i))
            if (withError=='kErrors'):
                #
                # Waiting for implemenation
                #
                print("waiting for implemeation:_variances")
                #reco.SetBinError (j, math.sqrt (abs(self._variances(i))))
            elif (withError=='kCovariance'):
                raise Exception('Error type not implemented')
                #reco.SetBinError (j, math.sqrt (abs (_cov(i,i))))
            elif (withError=='kCovToy'):
                raise Exception('Error type not implemented')
                #reco.SetBinError (j, math.sqrt (abs (_err_mat(i,i))));
        return reco
    
    def UnfoldWithErrors(self,withError,getWeights=False):

        if (not self._unfolded): 
            if (self._fail):
                return False
            rmeas= self._res.Hmeasured().Clone()

            if (self._meas.GetDimension() != rmeas.GetDimension() or
                self._meas.GetNbinsX()    != rmeas.GetNbinsX()    or
                self._meas.GetNbinsY()    != rmeas.GetNbinsY()    or
                self._meas.GetNbinsZ()    != rmeas.GetNbinsZ()):
                    raise Exception('data dimension different from response matrix reco dimension')
            self.Unfold()
            if (not self._unfolded):
                self._fail= True
                return False

        ok=True
        self._withError= withError

        if (getWeights and (withError=='kErrors' or withError=='kCovariance')):
            raise Exception('GetWeight not implemented') #shouldn't enter this case by current settings
            #if (not self._haveWgt):      
            #    GetWgt()
            #ok= self._haveWgt
        else: 
            if withError=='kErrors':
                if (not self._haveErrors):   
                    #
                    #waiting for implementation
                    #
                    print('Waiting for implementation: GetErrors()')
                    #self.GetErrors()
                ok=True
                #ok= self._haveErrors
            else:
                raise Exception('Other error type not implemented')

        if (not ok):
            self._fail= True
        return ok
    
    def Unfold(self):
        pass

    def Vmeasured(self):
    # Measured distribution as a vector.
        if (not self._vMes):
            self._vMes= RooUnfoldResponse.H2V(self._meas, self._res.GetNbinsMeasured(), self._overflow)
        return self._vMes



class RooUnfoldBayes(RooUnfold):
    def __init__(self, res, meas, niter):
        super(RooUnfoldBayes,self).__init__(res,meas)
        self._niter = niter
        self.Init_Bayes()

    def Init_Bayes(self):
        self._nc= self._ne= 0
        self._nbartrue= self._N0C= 0.0
        self.GetSettings_Bayes()
    
    def GetSettings_Bayes(self):
        self._minparm=1
        self._maxparm=15
        self._stepsizeparm=1
        self._defaultparm=4
    
    def Unfold(self):
        self.setup() #lower case version
        if (self._verbose >= 2):
            print('')
            #Print()
            #PrintMatrix(_Nji,"RooUnfoldBayes response matrix (Nji)")
        if (self._verbose >= 1):
            print("Now unfolding...")
        self.unfold()
        if (self._verbose >= 2):
            print('')

        self._rec = ROOT.TVectorD()
        self._rec.ResizeTo(self._nc)
        self._rec = self._nbarCi.Clone()
        self._rec.ResizeTo(self._nt) # drop fakes in final bin
        self._unfolded= True
        self._haveCov=  False

    @staticmethod
    def H2M (h, m, overflow):
    # TH2 -> TMatrixD
        if (not h): 
            return m
        first= 0 if overflow else 1
        nm= m.GetNrows()
        nt= m.GetNcols()
        for j in range (0,nt):
            for i in range(0,nm):
                m[i,j]= h.GetBinContent(i+first,j+first)
        return m


    def setup(self):
        self._nc = self._nt
        self._ne = self._nm
        self._nEstj = ROOT.TVectorD()
        self._nEstj.ResizeTo(self._ne)
        self._nEstj= self.Vmeasured().Clone()

        self._nCi = ROOT.TVectorD()
        self._nCi.ResizeTo(self._nt)
        self._nCi= self._res.Vtruth().Clone()

        self._Nji = ROOT.TMatrixD()
        self._Nji.ResizeTo(self._ne,self._nt)
        RooUnfoldBayes.H2M (self._res.Hresponse().Clone(), self._Nji, self._overflow)   # don't normalise, which is what _res->Mresponse() would give us

        if (self._res.FakeEntries()):
            fakes= self._res.Vfakes().Clone()
            nfakes= fakes.Sum()
            if (self._verbose>=0): 
                print("Add truth bin for ",nfakes, " fakes")
            self._nc+=1
            self._nCi.ResizeTo(self._nc)
            self._nCi[self._nc-1]= nfakes
            self._Nji.ResizeTo(self._ne,self._nc)
    
            for i in range (0, self._nm):
                self._Nji[i,self._nc-1]= fakes[i]
  
        self._nbarCi = ROOT.TVectorD()
        self._efficiencyCi = ROOT.TVectorD()
        self._Mij= ROOT.TMatrixD()
        self._P0C= ROOT.TVectorD()
        self._UjInv= ROOT.TVectorD()

        self._nbarCi.ResizeTo(self._nc)
        self._efficiencyCi.ResizeTo(self._nc)
        self._Mij.ResizeTo(self._nc,self._ne)
        self._P0C.ResizeTo(self._nc)
        self._UjInv.ResizeTo(self._ne)

        self._dnCidnEj = ROOT.TMatrixD()
        self._dnCidPjk = ROOT.TMatrixD()

        #ifndef OLDERRS
        if (self._dosys!=2):
            self._dnCidnEj.ResizeTo(self._nc,self._ne)
        #endif
        if (self._dosys):    
            self._dnCidPjk.ResizeTo(self._nc,self._ne*self._nc)

        #Initial distribution
        self._N0C= self._nCi.Sum()
        if (self._N0C!=0.0):
            self._P0C= self._nCi.Clone()
            self._P0C *= 1.0/self._N0C

    def getChi2(self,prob1, prob2, nevents):
    #calculate the chi^2. prob1 and prob2 are the probabilities
    #and nevents is the number of events used to calculate the probabilities
        chi2= 0.0
        n= prob1.GetNrows()
        if (self._verbose>=2):
            print("chi2 ",n," ",nevents) 
        for i in range(0,n):
            psum  = (prob1[i] + prob2[i])*nevents
            pdiff = (prob1[i] - prob2[i])*nevents
            if (psum > 1.0):
                chi2 = chi2 + (pdiff*pdiff)/psum
            else:
                chi2 = chi2 + (pdiff*pdiff)
        return chi2

    def unfold(self):
    # Calculate the unfolding matrix.
    # _niter = number of iterations to perform (3 by default).
    # _smoothit = smooth the matrix in between iterations (default false).
        PEjCi = ROOT.TMatrixD(self._ne,self._nc)
        PEjCiEff = ROOT.TMatrixD(self._ne,self._nc)
        for i in range(0,self._nc):
            if (self._nCi[i] <= 0.0): 
                self._efficiencyCi[i] = 0.0
                continue

            eff = 0.0
            for j in range(0,self._ne):
                response = self._Nji(j,i) / self._nCi[i]
                PEjCi[j,i] = PEjCiEff[j,i] = response  #efficiency of detecting the cause Ci in Effect Ej
                eff += response
        
            self._efficiencyCi[i] = eff
            effinv = 1.0/eff if eff > 0.0 else 0.0   #reset PEjCiEff if eff=0
            for j in range(0,self._ne):
                PEjCiEff[j,i] = PEjCiEff(j,i)*effinv
  

        PbarCi= ROOT.TVectorD(self._nc)

        for kiter in range(0, self._niter):
            
            if (self._verbose>=1):
                print("Iteration : %s"%kiter)

        #pdate prior from previous iteration
            if (kiter>0):
                self._P0C = PbarCi.Clone()
                self._N0C = self._nbartrue
    

            for j in range(0,self._ne): 
                Uj = 0.0
                for i in range(0, self._nc):
                    Uj += PEjCi(j,i) * self._P0C[i]
                self._UjInv[j] =  1.0/Uj if Uj > 0.0 else 0.0
            

            #Unfolding matrix M
            self._nbartrue = 0.0
            for i in range(0,self._nc):
                nbarC = 0.0
                for j in range(0,self._ne):
                    Mij = self._UjInv[j] * PEjCiEff(j,i) * self._P0C[i]
                    self._Mij[i,j]= Mij
                    nbarC += Mij * self._nEstj[j]
            
                self._nbarCi[i] = nbarC
                self._nbartrue += nbarC  # best estimate of true number of events
            

            # new estimate of true distribution
            PbarCi= self._nbarCi.Clone()
            PbarCi *= 1.0/self._nbartrue

        #ifndef OLDERRS
            if (self._dosys!=2): 
                if (kiter <= 0): 
                    self._dnCidnEj= self._Mij.Clone()
                else:
        #ifndef OLDMULT
                    en = ROOT.TVectorD(self._nc)
                    nr = ROOT.TVectorD(self._nc)
                    for i in range(0, self._nc):
                        #print("Check 0 nc0:",i,":",nr[i])
                        if (self._P0C[i]<=0.0): 
                            continue
                        ni= 1.0/(self._N0C*self._P0C[i])
                        en[i]= -ni*self._efficiencyCi[i]
                        nr[i]=  ni*self._nbarCi[i]
                        #print("Check 0 nc0 after update:",i,":",nr[i])
                    
                    M1= self._dnCidnEj.Clone()
                    M1.NormByColumn(nr,"M")
                    M2 = ROOT.TMatrixD(ROOT.TMatrixD.kTransposed, self._Mij)
                    M2.NormByColumn(self._nEstj,"M")
                    M2.NormByRow(en,"M")
                    M3 = ROOT.TMatrixD(M2, ROOT.TMatrixD.kMult, self._dnCidnEj)
                    self._dnCidnEj.Mult (self._Mij, M3)
                    self._dnCidnEj += self._Mij
                    self._dnCidnEj += M1
        #endif
        #endif

            #no need to smooth the last iteraction
            #if (_smoothit && kiter < (_niter-1)) smooth(PbarCi);

            #Chi2 based on Poisson errors
            chi2 = self.getChi2(PbarCi, self._P0C, self._nbartrue)
            if (self._verbose>=1):
                print("Chi^2 of change %s"%chi2)

            # and repeat
       



  

  
        

  

