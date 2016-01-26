direc=/data/jbrandstetter/CMGTools/rootFiles_160107/additionalSelection/
#wwwdirec=~/www/MVA_validation/mt/160107/
wwwdirec=~/www/MVA_validation/et/160107/


echo "Removing all existing .png files"
rm *.png
if [ ! -d "$wwwdirec" ];then 
    mkdir ${wwwdirec}
fi
if [ ! -d ${wwwdirec}"basis/" ];then 
    mkdir ${wwwdirec}"basis/"
fi
if [ ! -d ${wwwdirec}"basis_VBF/" ];then 
    mkdir ${wwwdirec}"basis_VBF/"
fi
if [ ! -d ${wwwdirec}"WJets_CR/" ];then 
    mkdir ${wwwdirec}"WJets_CR/"
fi



#python compare_bkg.py ${direc}"Ntuple_basis_MC_WJetsToLNu_madgraphMLM_mt.root" ${direc}"Ntuple_basis_QCD_datadriven_mt.root" ${direc}"Ntuple_basis_MC_TTbar_powheg_mt.root" ${direc}"Ntuple_basis_MC_DYJetsToLL_madgraphMLM_mt.root" ${direc}"Ntuple_basis_SingleMuon_Run2015D_mt.root" -t DYJets,TTbar,QCD,WJets,Data, -s basis

python compare_bkg.py ${direc}"Ntuple_basis_MC_WJetsToLNu_madgraphMLM_et.root" ${direc}"Ntuple_basis_QCD_datadriven_et.root" ${direc}"Ntuple_basis_MC_TTbar_powheg_et.root" ${direc}"Ntuple_basis_MC_DYJetsToLL_madgraphMLM_et.root" ${direc}"Ntuple_basis_SingleMuon_Run2015D_et.root" -t DYJets,TTbar,QCD,WJets,Data, -s basis

mv *.png ${wwwdirec}"basis/"

#python compare_bkg.py ${direc}"Ntuple_basis_VBF_MC_WJetsToLNu_madgraphMLM_mt.root" ${direc}"Ntuple_basis_VBF_QCD_datadriven_mt.root" ${direc}"Ntuple_basis_VBF_MC_TTbar_powheg_mt.root" ${direc}"Ntuple_basis_VBF_MC_DYJetsToLL_madgraphMLM_mt.root" ${direc}"Ntuple_basis_VBF_SingleMuon_Run2015D_mt.root" -t DYJets,TTbar,QCD,WJets,Data, -s basis_VBF

python compare_bkg.py ${direc}"Ntuple_basis_VBF_MC_WJetsToLNu_madgraphMLM_et.root" ${direc}"Ntuple_basis_VBF_QCD_datadriven_et.root" ${direc}"Ntuple_basis_VBF_MC_TTbar_powheg_et.root" ${direc}"Ntuple_basis_VBF_MC_DYJetsToLL_madgraphMLM_et.root" ${direc}"Ntuple_basis_VBF_SingleMuon_Run2015D_et.root" -t DYJets,TTbar,QCD,WJets,Data, -s basis_VBF

mv *.png ${wwwdirec}"basis_VBF/"

#python compare_bkg.py ${direc}"Ntuple_basis_mT70Cut_MC_WJetsToLNu_madgraphMLM_mt.root" ${direc}"Ntuple_basis_mT70Cut_QCD_datadriven_mt.root" ${direc}"Ntuple_basis_mT70Cut_MC_TTbar_powheg_mt.root" ${direc}"Ntuple_basis_mT70Cut_MC_DYJetsToLL_madgraphMLM_mt.root" ${direc}"Ntuple_basis_mT70Cut_SingleMuon_Run2015D_mt.root" -t DYJets,TTbar,QCD,WJets,Data, -s WJets_CR

python compare_bkg.py ${direc}"Ntuple_basis_mT70Cut_MC_WJetsToLNu_madgraphMLM_et.root" ${direc}"Ntuple_basis_mT70Cut_QCD_datadriven_et.root" ${direc}"Ntuple_basis_mT70Cut_MC_TTbar_powheg_et.root" ${direc}"Ntuple_basis_mT70Cut_MC_DYJetsToLL_madgraphMLM_et.root" ${direc}"Ntuple_basis_mT70Cut_SingleMuon_Run2015D_et.root" -t DYJets,TTbar,QCD,WJets,Data, -s WJets_CR

mv *.png ${wwwdirec}"WJets_CR/"


#python compare_bkg.py ${direc}"Ntuple_basis_VBF_MC_WJetsToLNu_madgraphMLM_mt_BDTscore.root" ${direc}"Ntuple_basis_VBF_QCD_datadriven_mt_BDTscore.root" ${direc}"Ntuple_basis_VBF_MC_TTbar_powheg_mt_BDTscore.root" ${direc}"Ntuple_basis_VBF_MC_DYJetsToLL_madgraphMLM_mt_BDTscore.root" ${direc}"Ntuple_basis_VBF_SingleMuon_Run2015D_mt_BDTscore.root" -t DYJets,TTbar,QCD,WJets,Data, -s basis_VBF

#mv *.png ${wwwdirec}"basis_VBF/"

#python compare_bkg.py ${direc}"Ntuple_basis_mT70Cut_MC_WJetsToLNu_madgraphMLM_mt_BDTscore.root" ${direc}"Ntuple_basis_mT70Cut_QCD_datadriven_mt_BDTscore.root" ${direc}"Ntuple_basis_mT70Cut_MC_TTbar_powheg_mt_BDTscore.root" ${direc}"Ntuple_basis_mT70Cut_MC_DYJetsToLL_madgraphMLM_mt_BDTscore.root" ${direc}"Ntuple_basis_mT70Cut_SingleMuon_Run2015D_mt_BDTscore.root" -t DYJets,TTbar,QCD,WJets,Data, -s WJets_CR

#mv *.png ${wwwdirec}"WJets_CR/"