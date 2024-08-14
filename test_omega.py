import sys
from openeye import oechem
from openeye import oeomega


def generate_conformers(smiles_list, output_file):
    omegaOpts = oeomega.OEOmegaOptions()

    print(oeomega.OEOmegaGetRelease())
    print(oeomega.OEOmegaGetLicensee())
    print(omegaOpts.GetSearchForceField())

    omegaOpts.SetRMSThreshold(0.2)
    omegaOpts.SetEnergyWindow(10)
    omegaOpts.SetMaxConfs(10000)
    omegaOpts.SetMaxSearchTime(600)

    omega = oeomega.OEOmega(omegaOpts)
    ofs = oechem.oemolostream()
    if not ofs.open(output_file):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % output_file)

    for smiles in smiles_list:
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smiles)
        ret_code = omega.Build(mol)
        if ret_code == oeomega.OEOmegaReturnCode_Success:
            conf_num = 0
            for conf in mol.GetConfs():
                oechem.OEWriteMolecule(ofs, conf)
                conf_num += 1
            print("Conformation number: "+str(conf_num))
        else:
            oechem.OEThrow.Warning("%s: %s" % (smiles, oeomega.OEGetOmegaError(ret_code)))


if __name__ == "__main__":
    smiles_list = ["O=C1NC(=O)c2c(NCC(=O)N3CCN(c4ccccc4)CC3)cccc21"]
    generate_conformers(smiles_list, output_file="test.pdb")
