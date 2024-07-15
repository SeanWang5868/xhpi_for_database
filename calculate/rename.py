def pi_rename(residue):
    if residue.name == 'TRP': 
        pi_residue_name = 'TRP_B'
    elif residue.name == 'HIS': 
        pi_residue_name = 'HIS'
    elif residue.name == 'TYR': 
        pi_residue_name = 'TYR'
    elif residue.name == 'PHE': 
        pi_residue_name = 'PHE'
    
    return pi_residue_name