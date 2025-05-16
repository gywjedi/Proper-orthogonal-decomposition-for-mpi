path = r"./"
for job_id in range(40,60):
    file_path = path + "/q_batch." + f"{job_id:03d}" +".sh"  # Specify the desired file path
    fname_id = 5000000 * job_id
    fname = "restart." + f"{fname_id}"
    fid = job_id
    fname_einstein = "profile_einstein" + f"{fid:02d}" +".dat"
    fname_ptensor = "./ptensor/ptensor.txt." + f"{fid:03d}"
    with open(file_path, "w") as file:
        # Write text to the file
        file.write("#!/bin/bash\n")
        file.write("#SBATCH -A chem\n")
        file.write("#SBATCH -p high_mem_cd\n")
        file.write("#SBATCH -N 1\n")
        file.write("#SBATCH --ntasks-per-node 36\n")
        file.write("#SBATCH -J x16\n")
        file.write("#SBATCH --mem=0\n")
        file.write("#SBATCH -c 1\n")
        file.write("#SBATCH -t 24:00:00\n")
        file.write("#SBATCH -o ./output.txt\n")
        file.write("#SBATCH -e ./cades-error.txt \n")
        file.write("\n")        
        file.write("cd $SLURM_SUBMIT_DIR\n")
        file.write("date\n")
        file.write("module load gcc\n")
        file.write("module load fftw\n")
        file.write("\n")
        file.write("\n")
        file.write("pwd\n")
        file.write("srun -n36 -c1 --cpu-bind=cores /home/gywjedi90/lammps_2023-12-15/build_2025-03-27/lmp_cades < in.nvt -var fname " + f"{fname} " + "-var fname_einstein " + f"{fname_einstein} " + "-var fname_ptensor " + f"{fname_ptensor} " + " &\n")
        file.write("wait")
