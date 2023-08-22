
# it is assumed that this script is run 
# in the same folder it resides in
work_dir=${PWD}
r_script="Rscript ${work_dir}/temp_bayes_prism_kang_test.R "
num_samp=1000
ncores=20

in_dir=${work_dir}/../data/single_cell_data/cibersort_kang/
out_dir=${work_dir}/../results/single_cell_data/bp_kang/


file_id_train=mono-kang_0
file_id_test=mono-kang_0


run_id="bp_kang"
lsf_file=${out_dir}/temp_bp_final.lsf
bsub -R "rusage[mem=48GB]" -W 24:00 -n 20 -q "normal" -o ${lsf_file} -J ${run_id} ${r_script} ${in_dir} ${out_dir} ${file_id_train} ${file_id_test} ${num_samp} ${ncores}
