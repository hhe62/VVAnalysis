#variables="pt mass zpt leppt dphiz1z2 drz1z2"
variables="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass3j Mass34j Mass4j MassFull Mass0jFull Mass1jFull Mass2jFull Mass3jFull Mass34jFull Mass4jFull"
#variables="dphiz1z2"
for var in $variables;do
  echo $var
  ./Utilities/scripts/saveUnfolded.py -a ZZ4l2016 -s LooseLeptons -l 35.9 -f ZZ4l2016 -sf data/scaleFactorsZZ4l2016.root -ls 2016fWUnc_full -vr ${var} --test --makeTotals --plotResponse --plotDir /afs/cern.ch/user/h/hehe/www/FullvarList_20May2021/UnfoldZZ4l2016_reg_MatrixPlots_20220308Check/ #--noNorm
done

echo "Job done================"
