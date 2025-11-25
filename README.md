# protein_SIP_pipeline
This is a Nextflow pipline that runs [Sipros4](https://github.com/thepanlab/Sipros4) and [IsopacketModeler](https://github.com/stavis1/isopacketModeler) to analyze protein stable isotope probing (SIP) data. 

## Installation
It is strongly recommended that this pipeline be run on high performance computing (HPC) resources as both component tools are substantially more resource intensive than standard proteomics pipelines. The provided example `nextflow.config` file is set up to work with the SLURM scheduler. See [footnote 1](#requirements) for more information.

Installing protein_SIP_pipeline requires three components: a copy of this Git repo, a working Nextflow installation, and a container runtime such as Apptainer or Docker.

### The protein_SIP_pipeline git repo:
To download the git repo either:
1. run `git clone https://github.com/stavis1/protein_SIP_pipeline.git`

or

2. click on the green `<> code` menu and then click `Download ZIP` and uncompress the resulting .zip file in the desired location.

### Nextflow:
To install Nextflow using conda:
run `conda env create -n nextflow_env -c bioconda nextflow=25.10.0`
If you wish to install nextflow without using conda please follow the installation instructions provided with the [nextflow documentation](https://www.nextflow.io/docs/latest/install.html). The example analysis and usage instructions will assume that conda was used for installation if not replace `conda run -n nextflow_env nextflow` with just `nextflow` in the example commands. 

### The container runtime:
The HPC environment that you use is likely to come with a preinstalled container runtime. For academic clusters this is likely to be Apptainer/Singularity.  The provided `nextflow.config` file assumes that Apptainer is the installed runtime. If not, then simply replace `apptainer` with the name of the installed runtime. Please talk to your HPC system administrator if you do not know if or which runtime is installed. 

You will need to be able to pull images from dockerhub on your HPC. Assuming you have access to Apptainer, try running the command `apptainer pull hello.sif docker://hello-world:latest` then checking that `hello.sif` has been created to test this ability. If you do not have this ability then follow these instructions: 
1. Download the git repo to a machine that does have this ability.
2. Navigate to the `cache` directory.
3. Run the `build_apptainer.sh` script.
4. Upload the `.img` files that were created to the `cache` directory in the git repo on your HPC

Other container runtimes will require slightly different pull commands, please consult with your HPC system administrator if you're using a different container runtime. 

Installing and configuring your own container runtime is outside of the scope of these install instructions so, if this is necessary, please follow the instructions in the documentation of your desired runtime. 

## Usage instructions
This pipeline uses a command line interface. Three configuration files are required to run the pipeline: `nextflow.config`, `config.cfg` and `design.tsv`. Templates for each of these files can be found in the `configs` directory. You will also need a Thermo .raw file, a converted .mzML file for each sample and a .fasta file database of protein sequences to search for.

To run the pipeline run:
`conda run -n nextflow_env nextflow -C nextflow.config /path/to/protein_SIP_pipeline/main.nf --design design.tsv --results_dir results/ -resume -with-report results/report.html`

### `nextflow.config`
This file should be set up once per HPC installation. It controls how Nextflow interacts with the HPC job scheduler and container runtime. It also controls the resources allocated to each job. The following instructions assume that your HPC uses SLURM and Apptainer. If not, please consult the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) to determine how to modify this file for your use case. You will need to know several pieces of information about your HPC environment and account, please consult your HPC's documentation or your system administrator to find these.

1. set the `queueSize` variable to the maximum number of jobs you are allowed to queue at once, minus one for the Nextflow parent job
2. set the `account` variable to the SLURM resource account you use when submitting SLURM jobs. This is likely not your user account. This value must be wrapped in quotes. 
3. set the `clusterOptions` variable to include the standard information you use to submit SLURM jobs. Typically you will need to change `partition_name` and `quality_of_service`. Leave `-N 1` as is, none of the child jobs submitted by the pipeline are able to use more than one node. This value must be wrapped in quotes.
4. The `cpus`, `memory`, and `time` variables for each job have been set to work well for several datasets run on the HPC I have access to. I suggest that you leave them unchanged unless you notice jobs failing due to resource limitations. If you wish to try optimizing these values see [footnote 2](#optimizing-resources). 

### `design.tsv`
This file specifies sample level information needed by the pipeline. It is a tab separated table of values. This file must be named `design.tsv`

 - `sample_ID` A sample identifier. It should not include spaces and must be globally unique within the run.
 - `raw_file` The thermo .raw filename for the sample. Only thermo .raw files are currently supported by Sipros.
 - `label_elm` The element that was used for isotopic labeling. Only one element is supported. Options are C, H, N, O, P, and S. Deuterium labeling (H) is supported but not advised, see [footnote 3](#deuterium-labeling). Leave this blank for unlabeled control samples, see [footnote 4](#unlabeled-controls). 
 - `label_integer` The mass number of the label isotope, e.g. 13 for ^13^C.  Leave this blank for unlabeled control samples. 
 - `config` The `config.cfg` file name to use for the sample. This file's name is not constrained.
 - `sipros_reduce` By default Sipros runs 100 database searches per file, one for each atom percent label incorporation from background to 100%. To run only the Nth search job set this value to N, e.g. to run only every other search set this to 2. 
 - `fasta` The fasta database to use for searching. 

### `config.cfg`
This file controls the settings for both Sipros and IsopacketModeler. Please consult the relevant documentation for each tool for detailed instructions regarding parameter settings. You will need at least one of these files for each run. If you wish to use different settings for different samples in the same run then one file will need to be made for each combination of settings. The most common use case for this feature is to enable searching different isotopic labels in the same run. 

## Example dataset  
An example dataset is coming soon. 

## Interpreting results  
The results folder will contain relevant outputs from both Sipros and IsopacketModeler. The important output from IsopacketModeler is `peptides.tsv` which contains information on all isotopically enriched peptides that pass the IsopacketModeler's quality control filters. This file is described in detail [here](https://github.com/stavis1/isopacketModeler). There will be one folder per sample that contains the Sipros outputs. These are described [here](https://github.com/thepanlab/Sipros4/wiki). 

## Footnotes
### requirements
This pipeline requires a POSIX-compatible OS (Linux, freeBSD, MacOS, etc.) running on x86 hardware. I have only tested this pipeline on the [CADES](https://www.ornl.gov/content/cades) cluster at ORNL so I do not know the compatibility with other systems. In principle this should be widely usable, as both Nextflow and the containers all tools run in are designed to be portable. Nextflow has executors that work with most common cluster and cloud schedulers such as AWS batch, Kubernetes, PBS, etc. Using these systems should only require modifications to the `nextflow.config` file. The containers run under both Docker and Apptainer runtimes. 

### optimizing resources
I would suggest ensuring that the values are a clean factor of a node's available resources as this will maximize resource utilization when a node is fully occupied by concurrent jobs. `large` is the most critical label to optimize, as this is used by the Sipros search step which runs 100 jobs per file and is typically the most computationally intensive. 

### deuterium labeling
The problem is that deuterium induces a retention time shift. This means that heavily labeled isotopologues will not co-elute with their light counterparts. IsopacketModeler assumes that isotopologues co-elute so this effect will distort the observed isotopic packet. The effect will be most severe in the context of the background + enriched models. This may not be a serious concern if you have good reason to believe that the isotope incorporation will be unimodal and reasonably uniform, i.e. the binomial model. This would be the case for labeling of species isolates in liquid culture.  I, however, do not have the data necessary to determine the severity of this problem. 

### unlabeled controls
IsopacketModeler uses unlabeled control samples to train a neural network classifier that filters out unlabeled peptides in the labeled samples. I suggest using this feature, as it both provides false discovery rate control on the identification of labeled peptides and speeds up the computation when many unlabeled peptides are present. Ideally the experiment is run with unlabeled control samples. However, if you are analyzing a dataset that does not have this feature it may be possible to use similar samples run on the same instrument with the same settings instead. I do not have the data necessary to benchmark how well this would work but, at least in principle, it is likely to give reasonable results. 

Sipros does allow for unlabeled database searches, however the settings are quite different between these two modes and the output files are incompatible. Currently, therefore, even unlabeled data must be searched using the labeled strategy. The same label isotope should be searched for in the unlabeled data as in the labeled data. If you are doing multiple label isotopes in the same experiment then the labeled files should be split among them. 

