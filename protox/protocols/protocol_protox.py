################ VERÓNICA GAMO PAREJO ##############################

# General imports 
import os
from urllib.request import urlopen


# Specific imports
from pyworkflow.protocol.params import PointerParam, EnumParam
from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem.utils import *
from protox import Plugin
import pandas as pd
from pwchem.constants import RDKIT_DIC


class ProtChemProtox(EMProtocol):

    """Toxicity prediction of small ligands with Protox"""
    
    _label = 'toxicity prediction of small ligands'
    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.smiles_list=[]


    def _defineParams(self, form):

        form.addSection(label='Input')
                      
        form.addParam('inputSet', PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Molecule for Toxicity Prediction', allowsNull=False,
                      help='Select the set of small molecules containing just one molecule for toxicity prediction.')
        
        form.addParam('predType', EnumParam, label='Prediction models', default=0,
                      choices = ["all",
                                "hepatotoxicity", "neurotoxicity", "nephrotoxicity", "respiratory toxicity", "cardiotoxicity",
                                "carcinogenicity", "immunotoxicity", "mutagenicity", "cytotoxicity", "BBB-barrier",
                                "ecotoxicity", "clinical toxicity", "nutritional toxicity", "nuclear receptor AhR", "nuclear receptor AR",
                                "nuclear receptor AR-LBD", " nuclear receptor aromatase", "nuclear receptor ER ", "nuclear receptor ER-LBD",
                                "nuclear receptor PPAR-Gamma", "nrf2/ARE", "HSE", "MMP",
                                "phosphoprotein p53", "ATAD5", "thyroid hormone receptor alpha", "thyroid hormone receptor beta",
                                "TTR", "RYR", "GABAR", "NMDAR",
                                "AMPAR", "KAR", "AChE", "CAR",
                                "PXR", "NADHOX", "VGSC", "NIS",
                                "CYP1A2", "CYP2C19", "CYP2C9", "CYP2D6",
                                "CYP3A4", "CYP2E1", "acute toxicity"
                            ],
                      help='Type of prediction to perform.')
        
        
    def _insertAllSteps(self):
        self._insertFunctionStep('extractSmile')
        self._insertFunctionStep('submitJob')
        self._insertFunctionStep('createOutputStep')

    def extractSmile(self):
        for mol in self.inputSet.get():
            smi = self.getSMI(mol, 1)
            self.smiles_list.append(smi)
            
    def getSMI(self, mol, nt):

        fnSmall = os.path.abspath(mol.getFileName())
        fnRoot, ext = os.path.splitext(os.path.basename(fnSmall))

        if ext != '.smi':
            outDir = os.path.abspath(self._getExtraPath())
            fnOut = os.path.abspath(self._getExtraPath(fnRoot + '.smi'))
            args = ' -i "{}" -of smi -o {} --outputDir {} -nt {}'.format(fnSmall, fnOut, outDir, nt)
            Plugin.runScript(self, 'rdkit_IO.py', args, env=RDKIT_DIC, cwd=outDir)    
        return self.parseSMI(fnOut)
        
    def parseSMI(self, smiFile):
        smi = None
        with open(smiFile) as f:
            for line in f:
                smi = line.split()[0].strip()
                if not smi.lower() == 'smiles':
                    break
        return smi
    
    def submitJob(self):
        models = ["acute_tox tox_targets ALL_MODELS",
    "dili", "neuro", "nephro", "respi", "cardio",
    "carcino", "immuno", "mutagen", "cyto", "bbb",
    "eco", "clinical", "nutri", "nr_ahr", "nr_ar",
    "nr_ar_lbd", "nr_aromatase", "nr_er", "nr_er_lbd",
    "nr_ppar_gamma", "sr_are", "sr_hse", "sr_mmp",
    "sr_p53", "sr_atad5", "mie_thr_alpha", "mie_thr_beta",
    "mie_ttr", "mie_ryr", "mie_gabar", "mie_nmdar",
    "mie_ampar", "mie_kar", "mie_ache", "mie_car",
    "mie_pxr", "mie_nadhox", "mie_vgsc", "mie_nis",
    "CYP1A2", "CYP2C19", "CYP2C9", "CYP2D6",
    "CYP3A4", "CYP2E1", "acute_tox tox_targets"]
        
        model = models[self.predType.get()]
        output_file = "results.csv"
        outDir = os.path.abspath(self._getExtraPath(output_file))
        Plugin.runProtox(self.smiles_list, model, outDir)
        print (f"Results saved in: {outDir}")

    def createOutputStep(self):
        csv_path = self._getExtraPath("results.csv")

        scores_df = pd.read_csv(csv_path, sep='\t', usecols=lambda col: col != 'Unnamed: 0')
        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='outputSmallMolecules')
        for mol in self.inputSet.get():
            fnSmall = os.path.abspath(mol.getFileName())
            smi = self.getSMI(mol, 1).strip()
            rows = scores_df[scores_df['input'].str.strip() == smi]
            print(f"Filas para {smi}:")
            print(rows)
            if not rows.empty:
                cid = self.getCIDFromSmiles(smi)
                name = self.getMainNameFromCID(cid)
                moleculeName = name if name is not None else smi

                for _, row in rows.iterrows():
                    smallMolecule = SmallMolecule(smallMolFilename=os.path.relpath(fnSmall), molName=moleculeName)
                    smallMolecule.type_Toxicity = pwobj.String(row['type'])
                    smallMolecule.target = pwobj.String(row['Target'])
                    smallMolecule.prediction = pwobj.Float(row['Prediction'])
                    smallMolecule.probability = pwobj.Float(row['Probability'])
                    outputSmallMolecules.append(smallMolecule)
                    
            
        outputSmallMolecules.updateMolClass()
        self._defineOutputs(outputSmallMolecules=outputSmallMolecules)

    def getCIDFromSmiles(self, smi):
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/cids/TXT" % smi
        try:
            with urlopen(url) as response:
                cid = response.read().decode('utf-8').split()[0]
        except Exception as e:
            #print(f"Error fetching CID for SMILES '{smi}': {e}")
            cid = None
        return cid
     
    def getMainNameFromCID(self,cid):
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/synonyms/TXT".format(cid)
        try:
            with urlopen(url) as response:
                r = response.read().decode('utf-8')
                synonyms = r.strip().split('\n')
                
                if synonyms:
                    main_name = synonyms[0].strip()
                else:
                    main_name = None
                
        except Exception as e:
            #print(f'Error fetching main name for CID {cid}: {e}')
            main_name = None
        
        return main_name
    

       

        
    
        

        
        

    
    
