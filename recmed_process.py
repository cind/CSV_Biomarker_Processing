import pandas as pd
import os
from pathlib import Path

# Create processed directory if it doesn't exist
processed_dir = Path('processed')
processed_dir.mkdir(exist_ok=True)

# Validate input file exists
file_path = Path('datasets/RECCMEDS_07Sep2024.csv')
if not file_path.exists():
    raise FileNotFoundError(f"Input file not found: {file_path}")

# Load the CSV file with error handling
try:
    data = pd.read_csv(file_path)
except Exception as e:
    print(f"Error reading CSV file: {e}")
    raise

# Extract relevant columns: 'RID' and 'CMMED'
subset_data = data[['RID', 'CMMED', 'VISDATE']]
subset_data = subset_data[subset_data['CMMED'] != '-4']

# Define multi-tiered medication taxonomy
# General Class > Subclass > Sub-Subclass
medication_taxonomy = {
    'Cardiovascular': {
        'Lipid-Lowering': ['lescol xl', 'lipitor', 'atorvastatin', 'simvastatin', 'rosuvastatin', 'pravastatin', 'lovastatin', 
                           'fenofibrate', 'gemfibrozil', 'niacin', 'ezetimibe', 'crestor', 'zocor', 'vytorin','zetia'],
        'Anti-Thrombotic': ['aspirin', 'clopidogrel', 'warfarin', 'heparin', 'xarelto', 'eliquis', 
                            'dabigatran', 'ticagrelor', 'prasugrel', 'plavix', 'coumadin', 'lovenox'],
        'Blood Pressure': ['amlodipine', 'doxazosin', 'lisinopril', 'losartan', 'losartin', 'metoprolol', 
                           'atenolol', 'cardura', 'hydrochlorothiazide', 'valsartan', 'ramipril', 'enalapril', 
                           'furosemide', 'spironolactone', 'diltiazem', 'verapamil', 'nifedipine', 'carvedilol', 
                           'propranolol', 'clonidine', 'irbesartan', 'norvasc', 'diovan', 'toprol xl', 'isosordide', 'imdur', 
                           'lasix', 'vasotec', 'benicar'],
    },
    'Metabolic': {
        'Thyroid': ['levothyroxine', 'synthroid', 'liothyronine', 'cytomel', 'methimazole', 'tapazole','propylthiouracil','levoxyl']
    },
    'Diabetic' : {
        'Diabetes Oral': ['metformin', 'glipizide', 'glyburide', 'sitagliptin', 'januvia', 'empagliflozin', 'jardiance', 'pioglitazone', 'glucophage', 'amaryl', 'avandia'],
        'Diabetes Injectable': ['insulin', 'liraglutide', 'victoza', 'exenatide', 'dulaglutide'],
    },
    'Psychiatric': {
        'Antidepressants': {
            'SSRI': ['lexapro', 'escitalopram', 'prozac', 'fluoxetine', 'sertraline', 'zoloft', 
                     'citalopram', 'celexa', 'paroxetine', 'paxil', 'fluvoxamine'],
            'SNRI': ['venlafaxine', 'effexor', 'duloxetine', 'cymbalta', 'desvenlafaxine', 
                     'pristiq', 'levomilnacipran', 'fetzima'],
            'Tricyclic': ['amitriptyline', 'nortriptyline', 'imipramine', 'desipramine', 
                          'doxepin', 'clomipramine'],
            'NDRI': ['bupropion', 'wellbutrin'],
            'Other': ['mirtazapine', 'trazodone', 'buspirone', 'lithium', 'lamotrigine', 'valproic acid']
        },
        'Benzodiazepines': ['alprazolam', 'xanax', 'diazepam', 'valium', 'clonazepam', 'klonopin', 
                            'lorazepam', 'ativan', 'temazepam', 'oxazepam'],
        'Antipsychotics': ['risperidone', 'abilify', 'olanzapine', 'zyprexa', 'quetiapine', 
                           'seroquel', 'aripiprazole', 'haloperidol', 'clozapine', 'ziprasidone']
    },
    'Pain Management': {
        'NSAIDs': ['ibuprofen', 'naproxen', 'indomethacin', 'advil', 'aleve', 'celecoxib', 
                  'diclofenac', 'meloxicam', 'piroxicam'],
        'Opioids': ['hydrocodone', 'oxycodone', 'morphine', 'codeine', 'tramadol', 
                   'fentanyl', 'dilaudid', 'methadone', 'buprenorphine'],
        'Analgesics': ['acetaminophen', 'tylenol', 'aspirin', 'asa'],
        'Corticosteroids': ['prednisone', 'hydrocortisone', 'dexamethasone', 
                            'methylprednisolone', 'betamethasone', 'triamcinolone']
    },
    'Supplements': {
        'Vitamins/Minerals': ['multi vitamin', 'omega 3', 'folic acid', 'vit c','vit e','mvi','vitamin a', 'vitamin b', 
                              'vitamin c', 'vitamin d', 'vitamin e', 'vitamin k', 'multivitamin', 
                              'b-12', 'b-complex', 'calcium', 'magnesium', 'zinc', 
                              'b-100 complex', 'iron', 'potassium', 'selenium', 'centrum silver', 'ocuvite', 'preservision', 
                              'slow mag', 'folbee','omega-3','coenzyme q10','b12','biotin'],
        'Herbal Supplements': ['lecithin', 'phosphatidyl choline', 'cod liver oil', 'l-tyrosine', 
                                'acetyl l-carnitine', 'bilberry', 'flax oil', 'lutein', 
                                'saw palmetto', 'fish oil', 'garlic', 'charcoal', 'turmeric', 
                                'ginkgo biloba', 'echinacea', 'glucosamine', 'chondroitin','flaxseed']
    },
    'Infectious Diseases': {
        'Antibiotics': ['augmentin','ciprofloxacin', 'bactrim', 'amoxicillin', 'penicillin', 
                        'doxycycline', 'azithromycin', 'zithromax', 'metronidazole', 
                        'cefazolin', 'sulfamethoxazole', 'levofloxacin', 'cephalexin', 
                        'clindamycin', 'vancomycin','cipro']
    },
    'Neurological': {
        'Anticonvulsants': ['gabapentin', 'pregabalin', 'topiramate', 'valproic acid', 
                            'carbamazepine', 'phenytoin', 'levetiracetam', 'lamotrigine'],
        'Parkinson\'s': ['carbidopa','levodopa', 'sinemet', 'pramipexole', 'ropinirole', 
                         'selegiline', 'amantadine']
    },
    'Gastrointestinal': {
        'GERD/PPI': ['omeprazole', 'pantoprazole', 'esomeprazole', 'lansoprazole', 
                     'nexium', 'prilosec', 'protonix', 'prevacid', 'aciphex', 'dexilant'],
        'H2 Blockers': ['ranitidine', 'famotidine', 'cimetidine', 'zantac', 'pepcid', 'tagamet'],
        'Other GI': ['colace','metoclopramide', 'ondansetron', 'loperamide', 'dicyclomine', 
                    'simethicone', 'bisacodyl', 'docusate', 'senna', 'fibercon', 'metamucil', 'immodium', 'zelnorm','miralax','tums']
    },
    'Respiratory': {
        'Inhalers': ['qvar','albuterol', 'fluticasone', 'budesonide', 'formoterol', 
                    'symbicort', 'advair', 'tiotropium', 'spiriva'],
        'Other Respiratory': ['montelukast', 'singulair', 'theophylline', 'ipratropium', 'zafirlukast']
    },
    'Urological': {
        'General Urological': ['tamsulosin', 'finasteride', 'oxybutynin', 'solifenacin', 
                                'tolterodine', 'flomax', 'proscar', 'detrol', 'dutasteride','terazosin','vesicare','avodart']
    },
    'Bone Health': {
        'Osteoporosis': ['alendronate', 'risedronate', 'ibandronate', 'zoledronic acid', 
                         'denosumab', 'teriparatide', 'raloxifene', 'fosamax','celebrex','actonel']
    },
    'Neurological Disorders': {
        'AD & Dementia': ['donepezil', 'aricept', 'memantine', 'namenda', 
                                    'rivastigmine', 'exelon', 'galantamine', 'razadyne', 
                                    'ebixa', 'cognex']
    },
    'Ophthalmologic': {
        'Glaucoma': ['xalatan', 'betoptics', 'cosopt', 'trusopt', 'travatan', 'alphagan'],
        'Eye Supplements': ['ocuvite', 'preservision']
    },
    'Allergy/Immunologic': {
        'Antihistamines': ['zyrtec', 'loratadine', 'allegra', 'benadryl', 'alovert','claritin'],
        'Nasal Sprays': ['flonase', 'nasacort', 'astelin']
    },
    'Sleep/Sedation': {
        'Sleep Aids': ['ambien', 'lunesta', 'melatonin'],
        'Anti-Anxiety': ['clonazepam']
    }
}

# Function to classify medications using multi-tiered taxonomy
def classify_medication_v8(med_name):
    med_name_lower = med_name.lower() if isinstance(med_name, str) else ''
    taxonomy = {'General_Class': 'Other', 'Subclass': 'Other', 'Sub_Subclass': 'Other'}
    
    for general_class, subclasses in medication_taxonomy.items():
        if isinstance(subclasses, dict):
            for subclass, sub_items in subclasses.items():
                if isinstance(sub_items, dict):
                    for sub_subclass, meds in sub_items.items():
                        if any(med in med_name_lower for med in meds):
                            taxonomy = {
                                'General_Class': general_class,
                                'Subclass': subclass,
                                'Sub_Subclass': sub_subclass
                            }
                            return taxonomy
                else:
                    if any(med in med_name_lower for med in sub_items):
                        taxonomy = {
                            'General_Class': general_class,
                            'Subclass': subclass,
                            'Sub_Subclass': 'None'
                        }
                        return taxonomy
        else:
            # Handle cases where subclasses might not be dictionaries
            pass
    
    return taxonomy

# Apply the grouping function to the 'CMMED' column
subset_data = subset_data.dropna(subset=['CMMED'])
classification = subset_data['CMMED'].apply(classify_medication_v8)
subset_data[['General_Class', 'Subclass', 'Sub_Subclass']] = pd.DataFrame(classification.tolist(), index=subset_data.index)

# Create mappings for each taxonomy level, ensuring 'Other' is always 0
general_class_mapping = {'Other': 0}
general_class_mapping.update({cls: i + 1 for i, cls in enumerate([x for x in subset_data['General_Class'].unique() if x != 'Other'])})

subclass_mapping = {'Other': 0}
subclass_mapping.update({sub: i + 1 for i, sub in enumerate([x for x in subset_data['Subclass'].unique() if x != 'Other'])})

sub_subclass_mapping = {'Other': 0}
sub_subclass_mapping.update({sub_sub: i + 1 for i, sub_sub in enumerate([x for x in subset_data['Sub_Subclass'].unique() if x != 'Other'])})

# Apply the mappings to encode the taxonomy
subset_data['General_Class_Encoded'] = subset_data['General_Class'].map(general_class_mapping)
subset_data['Subclass_Encoded'] = subset_data['Subclass'].map(subclass_mapping)
subset_data['Sub_Subclass_Encoded'] = subset_data['Sub_Subclass'].map(sub_subclass_mapping)

# Save the final dataset with the encoded columns
output_path = 'processed/encoded_medicines.csv'
subset_data['SCANDATE'] = subset_data['VISDATE'].copy()

# Save files with error handling
try:
    subset_data[['RID', 'CMMED', 'General_Class_Encoded', 'Subclass_Encoded', 
                 'Sub_Subclass_Encoded', 'SCANDATE']].to_csv(output_path, index=False)
    print(f"Successfully saved encoded data to {output_path}")
except Exception as e:
    print(f"Error saving encoded data: {e}")

# Save the mappings to separate CSV files
mapping_path_general = 'processed/general_class_mapping.csv'
pd.DataFrame(list(general_class_mapping.items()), columns=['General_Class', 'Encoding']).to_csv(mapping_path_general, index=False)

mapping_path_subclass = 'processed/subclass_mapping.csv'
pd.DataFrame(list(subclass_mapping.items()), columns=['Subclass', 'Encoding']).to_csv(mapping_path_subclass, index=False)

mapping_path_sub_subclass = 'processed/sub_subclass_mapping.csv'
pd.DataFrame(list(sub_subclass_mapping.items()), columns=['Sub_Subclass', 'Encoding']).to_csv(mapping_path_sub_subclass, index=False)

# Output the paths for all encoded files
print(f"Encoded CSV Output: {output_path}")
print(f"General Class Mapping CSV Output: {mapping_path_general}")
print(f"Subclass Mapping CSV Output: {mapping_path_subclass}")
print(f"Sub-Subclass Mapping CSV Output: {mapping_path_sub_subclass}")

# Count and print classification statistics
total_rows = len(subset_data)
classified_rows = len(subset_data[subset_data['General_Class_Encoded'] != 0])
unclassified_rows = total_rows - classified_rows

print("\nClassification Statistics:")
print(f"Total medications: {total_rows}")
print(f"Classified medications: {classified_rows} ({(classified_rows/total_rows*100):.1f}%)")
print(f"Unclassified medications: {unclassified_rows} ({(unclassified_rows/total_rows*100):.1f}%)")

# Add analysis of top unclassified medications
unclassified_meds = subset_data[subset_data['General_Class_Encoded'] == 0]['CMMED'].value_counts()
print("\nTop 10 most frequent unclassified medications:")
for med, count in unclassified_meds.head(10).items():
    print(f"{med}: {count} occurrences")