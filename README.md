# AFCluster, but modular

A modular reimplementation of [AF-Cluster](https://github.com/HWaymentSteele/AF_Cluster) for plug-n-play functionality, incorporation into computational workflows, and HPC/HTC deployment through ColabFold.

## Features

- **Modular Design**: Clean separation of MSA generation, clustering, and structure prediction
- **Slurm & Apptainer Compatibility**: Integration for high-performance and high-throughput computing clusters
- **Modular MSA Clustering**: Additional clustering methods can be easily swapped in.

## Installation/ Setup at CHTC at UW-Madison

0. There's no need to build a Docker, we provide scripts here under `scripts/chtc` that use the already-existing ColabFold Docker.

1. Git clone this repository to your home directory (i.e., `/home/waymentsteel/`)

2. [Just do this once, and again if your staging dir gets wiped]: download AF2 weights to local `staging` directory by calling `condor_submit download_af2_weights.sub` to run `download_af2_weights.sh`. Change username from waymentsteel to your own name.

3. Modify the condor submission script at `condor_submit_template.sub` to put in your username (look out for in paths too) and batch_name for your runs. This will launch jobs using the executable `run_afcluster_chtc.sh`. 

4. Create a file called `input_seqs.csv` that contains the IDs of the proteins you want to run in the first column, and the sequences in the second, with no header. for instance:

```
protein1,MKQVHHHSKLM
protein2,MAAALSDFERG
```

To submit a job, run  `condor_submit_template.sub`. This reads in the sequences in `input_seqs.csv`. You need to have the `condor_submit_template.sub`, `run_afcluster_chtc.sh`, and `input_seqs.csv` in the same dir. Results will also get stored there.

## Installation on other systems

### Prerequisites

- Python 3.8+
- ColabFold (for structure prediction)
- MMseqs2 (optional, for local MSA generation)

### Setup

1. **Clone the repository:**
   ```bash
   git clone git@github.com:gelnesr/AFCluster-2.git
   cd AFCluster-2

2. **Set up ColabFold locally:**
   
   If you do not already have ColabFold installed, please follow the installation instructions from [localcolabfold](https://github.com/YoshitakaMo/localcolabfold):
   
   ```bash
   # For Linux
   wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
   bash install_colabbatch_linux.sh
   ```
   After setting up ColabFold, set the relative path in the `configs/afcluster.yml` file. We recommend also setting up a cache directory.

3. **For HPC deployment, run the setup script:**

    This will automatically set up the enviornment and set the path for ColabFold to your `$SCRATCH/tools` folder. Edit this line in the .sh script to the appropriate directory. 
    ```bash
   bash scripts/env_setup.sh
   ```
4. **For Apptainer deployment, run the following set of commands.**
   
   First set up a .sif file to your `$SCRATCH/containers` folder and a cache directory at `$SCRATCH/cache`. Edit this line in the .sh script to the appropriate directory.
   
   ```bash
   bash scripts/build_apptainer.sh
   ```

    Then, run this which will automatically set up the enviornment and set the path for ColabFold to your `$SCRATCH/tools` folder. Edit this line in the .sh script to the appropriate directory. 

    ```bash
   bash scripts/env_setup.sh
   ```

   Then run `module load apptainer` or `module load singularity` to initialize an apptainer. You should set the `INPUT.fasta` in the command below before running:
   
   ```bash
   apptainer exec --nv \
        --bind "AFCluster-2:/w","$SCRATCH:$SCRATCH" \
        --env XDG_CACHE_HOME="$CACHE" \
        --env MPLCONFIGDIR="$CACHE" \
        "$IMG" bash -lc '
        cd /w
        source afc/bin/activate
        python afcluster.py --input INPUT.fasta
        '
    ```

    ****

## Usage

### Basic Usage

Run the pipeline on a FASTA file:

```bash
python afcluster.py --input sequences.fasta
```

## Parameters

- `--input`: Input FASTA file (required)
- `--msa`: Pre-computed MSA file (optional)
- `--config`: Config yaml file (see examples in `configs`)

### Configuration

Modify `configs/afcluster.yml` to adjust clustering parameters:

```yaml
keyword: "MAIN"           
gap_cutoff: 0.25 # or 0.1 for NeighborCluster protocol
random_seed: 42

cluster_method: "dbscan" # or "neighbor" for NeighborCluster
amber_relax: False

dbscan:
  min_samples: 10
  eps_val: null           
  min_eps: 3
  max_eps: 20.0
  eps_step: 0.5 

path_vars:
  PATH: "/path/to/colabfold/bin:$PATH"
  XDG_CACHE_HOME: "/path/to/cache"
  MPLCONFIGDIR: "/path/to/cache"
```

## Workflow

1. **Input Processing**: Load FASTA sequences
2. **MSA Generation**: Create multiple sequence alignments using `colabfold_batch` or local MMSeqs
3. **Sequence Filtering**: Remove sequences with high gap content
4. **Clustering**: Group similar sequences using specified method
5. **Structure Prediction**: Run ColabFold on each cluster
6. **Output**: Generate cluster-specific A3M files and predicted structures with corresponding json/png files

## Output Structure

```
output/
├── sequence_id/
│   ├── sequence_id.a3m          # Generated MSA
│   ├── clusters/
│   │   ├── sequence_id_000.a3m  # Cluster 0
│   │   ├── sequence_id_001.a3m  # Cluster 1
│   │   └── ...
│   └── preds/
│       ├── sequence_id_000/
│       │   ├── s0/              # Structure prediction seed 0
│       │   ├── s1/              # Structure prediction seed 1
│       │   └── ...
│       └── ...
```

## License

This project is based on the original [AF-Cluster](https://github.com/HWaymentSteele/AF_Cluster) implementation. 

## Citation

If you use AFCluster-2 in your research, please cite the following works and acknowledge this implementation:

```bibtex
@article{AFCluster,
  title={Predicting multiple conformations via sequence clustering and AlphaFold2},
  DOI={10.1038/s41586-023-06832-9},
  journal={Nature},
  author={Wayment-Steele, Hannah K. and Ojoawo, Adedolapo and Otten, Renee and Apitz, Julia M. and Pitsawong, Warintra and Hömberger, Marc and Ovchinnikov, Sergey and Colwell, Lucy and Kern, Dorothee},
  year={2023},
}
```

### Dependencies

This software builds upon the following tools and methods:

**AF-Cluster** - Wayment-Steele, H.K., Ojoawo, A., Otten, R. et al. Predicting multiple conformations via sequence clustering and AlphaFold2. *Nature* 625, 832–839 (2024). https://doi.org/10.1038/s41586-023-06832-9  
[[GitHub](https://github.com/HWaymentSteele/AF_Cluster)]

**ColabFold** - Mirdita, M., Schütze, K., Moriwaki, Y. et al. ColabFold: making protein folding accessible to all. *Nat Methods* 19, 679–682 (2022). https://doi.org/10.1038/s41592-022-01488-1  
[[GitHub](https://github.com/sokrypton/ColabFold)] | [[Local Installation](https://github.com/YoshitakaMo/localcolabfold)]

**MMseqs2** - Steinegger, M. & Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat Biotechnol* 35, 1026–1028 (2017). https://www.nature.com/articles/nbt.3988
[[GitHub](https://github.com/soedinglab/MMseqs2)]

**AlphaFold** - Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2
