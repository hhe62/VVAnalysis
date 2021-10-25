#variables="pt mass zpt leppt dphiz1z2 drz1z2"
#variables="dphiz1z2"
variables="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass3j Mass4j"
for var in $variables;do
  echo $var
  ./Utilities/scripts/plotUnfolded.py -a ZZ4l2017 -s LooseLeptons -l 41.5 -f ZZ4l2017 -vr ${var} --test --makeTotals --scaleymin 0.3 --scaleymax 1.2 --unfoldDir /afs/cern.ch/user/h/hehe/www/FullvarList_20May2021/UnfoldZZ4l2017/
  #--plotResponse
done