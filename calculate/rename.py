def pi_rename(residue,ring_length):

    if residue.name == 'TRP' and ring_length == 6: 
        pi_residue_name = 'TRP_B'
    elif residue.name == 'HIS': 
        pi_residue_name = 'HIS'
    elif residue.name == 'TYR': 
        pi_residue_name = 'TYR'
    elif residue.name == 'PHE': 
        pi_residue_name = 'PHE'
    else:
        pi_residue_name = 'TRP_A'
    
    return pi_residue_name