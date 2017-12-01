
# coding: utf-8

# In[3]:


# # # Requirements # # # 

import numpy as np 

# # # Constants # # # # 

# Import peptides and TOP from R script
# dict_peptide_positive = {'patient': 'peptides'}
# dict_peptide_negative =  {'patient': 'peptides'}
# TOP = top List For one allele And one k-fold 

# # # Functions # # # 

# Classification by occurence of allele-specific peptides (for each patient & threshold n)
def classify(allele, n, sequences,TOP): 
            counts = 0
            for peptide in TOP:
                if peptide in sequences:
                    counts += 1
            if counts >= n:
                return True
            if counts < n:
                return False
           
            
# Evaluation of output dependent on patient (allele-positive/negative)
def evaluate_p(patient, Allele_Patients_positive, output):
    TP, FN = 0, 0
    for patient in Allele_Patients_positive: 
        if output:
            TP = 1 
        if not output:
            FN = 1
    return TP, FN

def evaluate_n(patient, Allele_Patients_negative, output):
    FP, TN = 0, 0
    for patient in Allele_Patients_negative:
        if output:
            FP = 1
        if not output:
            TN = 1
    return FP, TN

def ROC (x, y): 
    xx = np.asarray(x)
    yy = np.asarray(y)
    AUC = np.trapz(yy, xx)
    return AUC

# # # Main # # # 

def main(TOP, dict_peptide_positive, dict_peptide_negative, allele, n):            # Define all thresholds to be parsed
 
    Allele_Patients_positive = dict_peptide_positive.keys()
    Allele_Patients_negative = dict_peptide_negative.keys()

    TP, FN, FP, TN = 0, 0, 0, 0
    for patient in dict_peptide_positive.keys():
        sequences = dict_peptide_positive[patient] 
        output = classify(allele,n,sequences, TOP)
        TP_n, FN_n = evaluate_p(patient, Allele_Patients_positive, output)
        TP += TP_n
        FN += FN_n
    for patient in dict_peptide_negative.keys():
        sequences = dict_peptide_negative[patient] 
        output = classify(allele,n,sequences, TOP)
        FP_n, TN_n = evaluate_n(patient, Allele_Patients_negative, output)
        TN += TN_n
        FP += FP_n
        
    return TP,FN,TN,FP

