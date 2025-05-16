from rdkit import Chem
from rdkit.Chem import Descriptors
from copy import deepcopy

# ここは実際のAcceptors, AliphaticRings, StructuralAlertsに置き換えてください
Acceptors = []  # 例: [Chem.MolFromSmarts('...'), ...]
AliphaticRings = Chem.MolFromSmarts('[R0]')  # 適切なSMARTSにしてください
StructuralAlerts = []  # 例: [Chem.MolFromSmarts('...'), ...]

def qed(mol):
    matches = []
    if mol is None:
        raise ValueError("mol argument is None")
    x = [0] * 9
    x[0] = Descriptors.MolWt(mol)
    x[1] = Descriptors.MolLogP(mol)
    for hba in Acceptors:
        if mol.HasSubstructMatch(hba):
            matches = mol.GetSubstructMatches(hba)
            x[2] += len(matches)
    x[3] = Descriptors.NumHDonors(mol)          
    x[4] = Descriptors.TPSA(mol)
    x[5] = Descriptors.NumRotatableBonds(mol)
    x[6] = len(Chem.GetSSSR(Chem.DeleteSubstructs(deepcopy(mol), AliphaticRings)))#lenで環の長さを取得
    for alert in StructuralAlerts:
        if mol.HasSubstructMatch(alert):
            x[7] += 1
    ro5_failed = 0
    if x[3] > 5:
        ro5_failed += 1
    if x[2] > 10:
        ro5_failed += 1
    if x[0] >= 500:
        ro5_failed += 1
    if x[1] > 5:
        ro5_failed += 1
    x[8] = ro5_failed

    for i, xi in enumerate(x):
        print(f"x[{i}]: value={xi}, type={type(xi)}")

if __name__ == "__main__":
    test_smiles = "CCO"  # ここを好きなSMILESに変えてテスト可能
    mol = Chem.MolFromSmiles(test_smiles)
    qed(mol)
