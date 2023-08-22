
# it is assumed that this script is run 
# in the same folder it resides in
work_dir=${PWD}
r_script="Rscript ${work_dir}/bayes_prism_liver_test.R "
num_samp=1000
ncores=20

in_dir=${work_dir}/../data/single_cell_data/cibersort_liver/
out_dir=${work_dir}/../results/single_cell_data/bp_liver/


file_id_train=all-liver_99
file_id_test=all-liver_99


run_id="bp_liver"
lsf_file=${out_dir}/temp_bp_final.lsf
bsub -R "rusage[mem=48GB]" -W 24:00 -n 20 -q "normal" -o ${lsf_file} -J ${run_id} ${r_script} ${in_dir} ${out_dir} ${file_id_train} ${file_id_test} ${num_samp} ${ncores}
