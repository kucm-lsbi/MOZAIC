import os
import pandas as pd
import logging

def save_xlsx(bank, output_path):

    if not bank:
        logging.warning("No bank to save.")
        return

    data_rows = []
    
    for res in bank:
        scores = res.get('scores', {})
        
        row = {
            'SMILES': res.get('product_smiles'),
            'Affinity': scores.get('affinity'),
            'QED': scores.get('qed'),
            'SA': scores.get('sa'),
            'Objective': scores.get('objective', 0.0)
        }
        
        # Reaction History (Step 1, Step 2 ...)
        rxn_history = res.get('rxn_history', [])
        for i, step in enumerate(rxn_history):
            prefix = f"Step {i+1}"
            row[f"{prefix} FG"]       = step.get('functional_group')
            row[f"{prefix} Position"] = step.get('rxn_position')
            row[f"{prefix} Fragment"] = step.get('selected_fragment')
            row[f"{prefix} Product"]  = step.get('reacted_compound')
            
        data_rows.append(row)

    df = pd.DataFrame(data_rows)

    # Align
    if 'Objective' in df.columns:
        df = df.sort_values(by='Objective', ascending=False)

    numeric_cols = ['Affinity', 'QED', 'SA', 'Objective']
    avg_data = {'SMILES': 'AVERAGE'}
    
    for col in numeric_cols:
        if col in df.columns:
            avg_val = df[col].mean()
            avg_data[col] = avg_val
            
    avg_df = pd.DataFrame([avg_data])
    df = pd.concat([df, avg_df], ignore_index=True)

    df = df.round(3)

    base_columns = ['SMILES', 'Objective', 'Affinity', 'QED', 'SA']
    
    step_columns = [c for c in df.columns if c not in base_columns]
    
    # Step 1 -> Step 2 ...
    final_columns = base_columns + step_columns
    final_columns = [c for c in final_columns if c in df.columns]
    
    df = df[final_columns]

    # To excel
    try:
        df.to_excel(output_path, index=False)
        
        avg_aff = avg_data.get('Affinity', 0)
        avg_qed = avg_data.get('QED', 0)
        avg_sa  = avg_data.get('SA', 0)
        
        logging.info(f"Results saved to {output_path}\n")
        logging.info(f"Average Stats -> Affinity: {avg_aff:.3f}, QED: {avg_qed:.3f}, SA: {avg_sa:.3f}\n")
        
    except Exception as e:
        logging.error(f"Failed to save Excel file: {e}")

def save_pose(bank, pose_dir, receptor_path):

    # rank_1.pdb, rank_2.pdb ... (complex pdb)

    if not bank:
        logging.warning("No poses to save.")
        return

    if not os.path.exists(pose_dir):
        os.makedirs(pose_dir)

    # Receptor 
    receptor_lines = []
    try:
        with open(receptor_path, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM', 'CONECT')):
                    receptor_lines.append(line)
    except Exception as e:
        logging.error(f"Failed to read receptor file: {e}")
        return

    sorted_bank = sorted(
        bank, 
        key=lambda x: x.get('scores', {}).get('objective', 0.0), 
        reverse=True
    )

    # Ranking
    count = 0
    for idx, res in enumerate(sorted_bank, start=1):
        # Ligand pose PDBQT
        pose_str = res.get('scores', {}).get('pose')
        
        if not pose_str:
            continue

        output_filename = os.path.join(pose_dir, f"rank_{idx}.pdb")
        
        try:
            # Receptor -> TER -> Ligand -> END
            with open(output_filename, 'w') as f:
                f.writelines(receptor_lines)
                
                f.write("TER\n")
                
                for line in pose_str.splitlines():
                    if line.startswith(('ATOM', 'HETATM')):
                        f.write(line + "\n")

                f.write("END\n")
                
            count += 1
            
        except Exception as e:
            logging.warning(f"Failed to save pose rank {idx}: {e}")

    logging.info(f"Saved {count} complex PDB files to {pose_dir}\n")