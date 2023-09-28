{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "ad233370",
   "metadata": {},
   "source": [
    "# scVelo Analysis on Myo combined dataset\n",
    "## Author: Lukas Tombor\n",
    "## Date: 22-02-17\n",
    "\n",
    "This analysis follows a tutorial based [here](https://scvelo.readthedocs.io/getting_started/)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "d1dd12ba",
   "metadata": {},
   "source": [
    "### 1- Load libraries"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "3eefb0e7",
   "metadata": {},
   "outputs": [],
   "source": [
    "import scvelo as scv\n",
    "import scanpy\n",
    "import os\n",
    "import pandas as pd\n",
    "import glob\n",
    "import numpy as np\n",
    "import scanpy as sc\n",
    "import scanorama\n",
    "import re\n",
    "\n",
    "\n",
    "scv.set_figure_params()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "0b170a64",
   "metadata": {},
   "source": [
    "### 2-Load datasets\n",
    "\n",
    "Use the output .loom files from velocyto which are located here: [scStorage/Lukas/Analysis/EHT/Velocity/out](files/../Velocity/out/)\n",
    "\n",
    "Script for running velocyto: [script](files/../Velocity/run_velocity_seq.sh)\n",
    "\n",
    "Commands: [Commands.txt](files/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Velocity/Commands.txt)\n",
    "          [rawcommands.txt](/../outcommands.txt)\n",
    "\n",
    "Velocyto Tutorial: [tutorial](http://velocyto.org/)\n",
    "\n",
    "\n",
    "For loop for all .loom files:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "f7ff7f85",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "[AnnData object with n_obs × n_vars = 6882 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 3527 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 1967 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 4014 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 1670 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 5671 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 1624 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 5928 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 1305 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 5924 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 4544 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 1007 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 1030 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 1786 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 2214 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 1744 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', AnnData object with n_obs × n_vars = 3077 × 28002\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced']\n"
     ]
    }
   ],
   "source": [
    "looms=list()\n",
    "names=list()\n",
    "mypath=\"/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Velocity/out/Young\"\n",
    "\n",
    "for f in glob.glob(glob.escape(mypath) + \"/*.loom\"):\n",
    "    adata=scv.read(f,cache=True)\n",
    "    adata.var_names_make_unique()\n",
    "    looms.append(adata)\n",
    "    names.append(f)\n",
    "print(looms)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "7615dd37",
   "metadata": {
    "scrolled": false
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "MF17010\n",
      "MF17013\n",
      "S64-104182-006\n",
      "MF17014\n",
      "S67-104182-007\n",
      "MF17015\n",
      "S68-104182-008\n",
      "MF17016\n",
      "S69-104182-009\n",
      "MF17017\n",
      "MF17018\n",
      "d1-103471-001-004\n",
      "d3-103471-001-002\n",
      "d7-103471-001-001\n",
      "d14-103548-001-005\n",
      "Hom-103471-001-003\n",
      "d28-103548-001-001\n"
     ]
    }
   ],
   "source": [
    "#getsamplenames \n",
    "samplenames=list()\n",
    "for n in names:\n",
    "    samplename = re.search(\"/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Velocity/out/Young/(.+?)_\", n).group(1)\n",
    "    print(samplename)\n",
    "    samplenames.append(samplename)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "adb8be26",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "                   Timepoint    Age   Batch\n",
      "MF17010                  Hom  Young  Batch1\n",
      "MF17013                   d1  Young  Batch1\n",
      "MF17014                   d3  Young  Batch1\n",
      "MF17015                   d5  Young  Batch1\n",
      "MF17016                   d7  Young  Batch1\n",
      "MF17017                  d14  Young  Batch1\n",
      "MF17018                  d28  Young  Batch1\n",
      "d7-103471-001-001         d7  Young  Batch2\n",
      "d3-103471-001-002         d3  Young  Batch2\n",
      "Hom-103471-001-003       Hom  Young  Batch2\n",
      "d1-103471-001-004         d1  Young  Batch2\n",
      "d28-103548-001-001       d28  Young  Batch3\n",
      "d14-103548-001-005       d14  Young  Batch3\n",
      "S64-104182-006           d28  Young  Batch4\n",
      "S67-104182-007           d14  Young  Batch4\n",
      "S68-104182-008            d7  Young  Batch4\n",
      "S69-104182-009           Hom  Young  Batch4\n"
     ]
    }
   ],
   "source": [
    "\n",
    "anno=pd.read_csv(\"/media/tombor/Helios_scStorage/Lukas/Analysis/EHT/Velocity/Anno.csv\", sep = \",\", index_col=0)\n",
    "print(anno)\n",
    "for i in range(len(looms)):\n",
    "    looms[i].obs[\"Sample\"] = samplenames[i]\n",
    "    looms[i].obs[\"Timepoint\"] = anno[\"Timepoint\"][samplenames[i]]\n",
    "    looms[i].obs[\"Age\"] = anno[\"Age\"][samplenames[i]]\n",
    "    looms[i].obs[\"Batch\"] = anno[\"Batch\"][samplenames[i]]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "8e91939c",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "                                        EGFP\n",
      "CellID                                      \n",
      "S69_104182-009_R9VNF:ACCACAATCACAATGCx   0.0\n",
      "S69_104182-009_R9VNF:AATTTCCGTCCAGTTAx   0.0\n",
      "S69_104182-009_R9VNF:AAAGTCCCACTAGAGGx   0.0\n",
      "S69_104182-009_R9VNF:ACCAACACACTTGTCCx   0.0\n",
      "S69_104182-009_R9VNF:ACCGTTCGTACTCAACx   0.0\n",
      "...                                      ...\n",
      "S69_104182-009_R9VNF:TTTGGAGGTCATCGGCx   0.0\n",
      "S69_104182-009_R9VNF:TTTGGAGCAACGTTACx   0.0\n",
      "S69_104182-009_R9VNF:TTTGTTGAGAGGTCACx   0.0\n",
      "S69_104182-009_R9VNF:TTTGTTGAGTAATTGGx   0.0\n",
      "S69_104182-009_R9VNF:TTTGGTTTCTCTCTTCx   0.0\n",
      "\n",
      "[1305 rows x 1 columns]\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "0.0"
      ]
     },
     "execution_count": 5,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "# Identify GFP pos cells\n",
    "\n",
    "df = sc.get.obs_df(looms[8], keys = [\"EGFP\"])\n",
    "print(df)\n",
    "\n",
    "s = pd.Series(df[\"EGFP\"])\n",
    "s.sum()\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "4b5b8ce9",
   "metadata": {},
   "source": [
    "Basic preprocessing"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "b6ab4578",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "(6882, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(3527, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(1967, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(4014, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(1670, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(5671, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(1624, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(5928, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(1305, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(5924, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(4544, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(1007, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(1030, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(1786, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(2214, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(1744, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n",
      "(3077, 28002)\n",
      "Filter cells\n",
      "MT-filter\n",
      "log1p\n",
      "HVG calc\n",
      "Regression\n",
      "Scale\n",
      "save to list...\n"
     ]
    }
   ],
   "source": [
    "processed_adatas=list()\n",
    "\n",
    "for adata in looms:\n",
    "    print(adata.shape)\n",
    "    print(\"Filter cells\")\n",
    "    sc.pp.filter_cells(adata, min_genes=200)\n",
    "    print(\"MT-filter\")\n",
    "    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'\n",
    "    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)\n",
    "    print(\"log1p\")\n",
    "    sc.pp.log1p(adata)\n",
    "    print(\"HVG calc\")\n",
    "    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\n",
    "    print(\"Regression\")\n",
    "    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'], n_jobs = 12)\n",
    "    print(\"Scale\")\n",
    "    sc.pp.scale(adata, max_value=10)\n",
    "    print(\"save to list...\")\n",
    "    processed_adatas.append(adata)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "a3c072fc",
   "metadata": {},
   "source": [
    "Combine adatas with scanorama "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "e223741b",
   "metadata": {
    "scrolled": false
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "AnnData object with n_obs × n_vars = 1962 × 28002\n",
      "    obs: 'Sample', 'Timepoint', 'Age', 'Batch', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'mean', 'std'\n",
      "    uns: 'log1p', 'hvg'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced'\n"
     ]
    }
   ],
   "source": [
    "print(processed_adatas[2])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "7399a0f1",
   "metadata": {
    "scrolled": false
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "476034\n",
      "28002\n"
     ]
    }
   ],
   "source": [
    "# Intersect highly variable genes\n",
    "\n",
    "def flatten(t):\n",
    "    return [item for sublist in t for item in sublist]\n",
    "\n",
    "hvar = list()\n",
    "for adata in processed_adatas:\n",
    "    hvar.append(adata.var.highly_variable.index.tolist())\n",
    "\n",
    "hvar = flatten(hvar)\n",
    "\n",
    "print(len(hvar))\n",
    "\n",
    "c = pd.Series(hvar).value_counts()\n",
    "\n",
    "\n",
    "# hvars in at least 3 datasets\n",
    "hvars = c[c > 2]\n",
    "hvars = list(hvars.index)\n",
    "\n",
    "\n",
    "print(len(hvars))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "4bdc9781",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Performing cosine normalization...\n",
      "Starting MNN correct iteration. Reference batch: 0\n",
      "Step 1 of 16: processing batch 1\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 2 of 16: processing batch 2\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 3 of 16: processing batch 3\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 4 of 16: processing batch 4\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 5 of 16: processing batch 5\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 6 of 16: processing batch 6\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 7 of 16: processing batch 7\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 8 of 16: processing batch 8\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 9 of 16: processing batch 9\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 10 of 16: processing batch 10\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 11 of 16: processing batch 11\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 12 of 16: processing batch 12\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 13 of 16: processing batch 13\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 14 of 16: processing batch 14\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 15 of 16: processing batch 15\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "Step 16 of 16: processing batch 16\n",
      "  Looking for MNNs...\n",
      "  Computing correction vectors...\n",
      "  Removing components...\n",
      "  Adjusting variance...\n",
      "  Applying correction...\n",
      "MNN correction complete. Gathering output...\n",
      "Packing AnnData object...\n",
      "Done.\n",
      "(AnnData object with n_obs × n_vars = 53699 × 28002\n",
      "    obs: 'Sample', 'Timepoint', 'Age', 'Batch', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand', 'mt', 'n_cells_by_counts-0', 'mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'highly_variable-0', 'means-0', 'dispersions-0', 'dispersions_norm-0', 'mean-0', 'std-0', 'n_cells_by_counts-1', 'mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'highly_variable-1', 'means-1', 'dispersions-1', 'dispersions_norm-1', 'mean-1', 'std-1', 'n_cells_by_counts-10', 'mean_counts-10', 'pct_dropout_by_counts-10', 'total_counts-10', 'highly_variable-10', 'means-10', 'dispersions-10', 'dispersions_norm-10', 'mean-10', 'std-10', 'n_cells_by_counts-11', 'mean_counts-11', 'pct_dropout_by_counts-11', 'total_counts-11', 'highly_variable-11', 'means-11', 'dispersions-11', 'dispersions_norm-11', 'mean-11', 'std-11', 'n_cells_by_counts-12', 'mean_counts-12', 'pct_dropout_by_counts-12', 'total_counts-12', 'highly_variable-12', 'means-12', 'dispersions-12', 'dispersions_norm-12', 'mean-12', 'std-12', 'n_cells_by_counts-13', 'mean_counts-13', 'pct_dropout_by_counts-13', 'total_counts-13', 'highly_variable-13', 'means-13', 'dispersions-13', 'dispersions_norm-13', 'mean-13', 'std-13', 'n_cells_by_counts-14', 'mean_counts-14', 'pct_dropout_by_counts-14', 'total_counts-14', 'highly_variable-14', 'means-14', 'dispersions-14', 'dispersions_norm-14', 'mean-14', 'std-14', 'n_cells_by_counts-15', 'mean_counts-15', 'pct_dropout_by_counts-15', 'total_counts-15', 'highly_variable-15', 'means-15', 'dispersions-15', 'dispersions_norm-15', 'mean-15', 'std-15', 'n_cells_by_counts-16', 'mean_counts-16', 'pct_dropout_by_counts-16', 'total_counts-16', 'highly_variable-16', 'means-16', 'dispersions-16', 'dispersions_norm-16', 'mean-16', 'std-16', 'n_cells_by_counts-2', 'mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'highly_variable-2', 'means-2', 'dispersions-2', 'dispersions_norm-2', 'mean-2', 'std-2', 'n_cells_by_counts-3', 'mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'highly_variable-3', 'means-3', 'dispersions-3', 'dispersions_norm-3', 'mean-3', 'std-3', 'n_cells_by_counts-4', 'mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'highly_variable-4', 'means-4', 'dispersions-4', 'dispersions_norm-4', 'mean-4', 'std-4', 'n_cells_by_counts-5', 'mean_counts-5', 'pct_dropout_by_counts-5', 'total_counts-5', 'highly_variable-5', 'means-5', 'dispersions-5', 'dispersions_norm-5', 'mean-5', 'std-5', 'n_cells_by_counts-6', 'mean_counts-6', 'pct_dropout_by_counts-6', 'total_counts-6', 'highly_variable-6', 'means-6', 'dispersions-6', 'dispersions_norm-6', 'mean-6', 'std-6', 'n_cells_by_counts-7', 'mean_counts-7', 'pct_dropout_by_counts-7', 'total_counts-7', 'highly_variable-7', 'means-7', 'dispersions-7', 'dispersions_norm-7', 'mean-7', 'std-7', 'n_cells_by_counts-8', 'mean_counts-8', 'pct_dropout_by_counts-8', 'total_counts-8', 'highly_variable-8', 'means-8', 'dispersions-8', 'dispersions_norm-8', 'mean-8', 'std-8', 'n_cells_by_counts-9', 'mean_counts-9', 'pct_dropout_by_counts-9', 'total_counts-9', 'highly_variable-9', 'means-9', 'dispersions-9', 'dispersions_norm-9', 'mean-9', 'std-9'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced', [0,        new cell  ref cell  original batch\n",
      "0             0      5846               0\n",
      "1             0      1537               0\n",
      "2             0      1258               0\n",
      "3             1      3606               0\n",
      "4             1      3116               0\n",
      "...         ...       ...             ...\n",
      "31953      3526      4227               0\n",
      "31954      3526      4705               0\n",
      "31955      3526      2827               0\n",
      "31956      3526       287               0\n",
      "31957      3526      4308               0\n",
      "\n",
      "[31958 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0      5418               0\n",
      "1             0      3310               0\n",
      "2             0      1491               0\n",
      "3             0      3029               0\n",
      "4             0       667               0\n",
      "...         ...       ...             ...\n",
      "21855      1960      3088               0\n",
      "21856      1960       304               0\n",
      "21857      1960      4431               0\n",
      "21858      1960      2171               0\n",
      "21859      1960      7618               1\n",
      "\n",
      "[21860 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0      8259               1\n",
      "1             1      8577               1\n",
      "2             1      6913               1\n",
      "3             1      7487               1\n",
      "4             1      7540               1\n",
      "...         ...       ...             ...\n",
      "25076      4011      9188               1\n",
      "25077      4011      5306               0\n",
      "25078      4011     11862               2\n",
      "25079      4011     10953               2\n",
      "25080      4011     12098               2\n",
      "\n",
      "[25081 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             1     12002               2\n",
      "1             1     10856               2\n",
      "2             1     10945               2\n",
      "3             1     10700               2\n",
      "4             1     10670               2\n",
      "...         ...       ...             ...\n",
      "20923      1667     10955               2\n",
      "20924      1667     11226               2\n",
      "20925      1667     12127               2\n",
      "20926      1667     12154               2\n",
      "20927      1667     11045               2\n",
      "\n",
      "[20928 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     17179               4\n",
      "1             0     11461               2\n",
      "2             1     14150               3\n",
      "3             1     12762               3\n",
      "4             1     17746               4\n",
      "...         ...       ...             ...\n",
      "30234      5667     13946               3\n",
      "30235      5667     12439               3\n",
      "30236      5667     15745               3\n",
      "30237      5667     12504               3\n",
      "30238      5667     14334               3\n",
      "\n",
      "[30239 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0      9383               1\n",
      "1             1     10956               2\n",
      "2             1     11739               2\n",
      "3             1     10542               2\n",
      "4             1     11353               2\n",
      "...         ...       ...             ...\n",
      "19106      1623     10608               2\n",
      "19107      1623     16781               4\n",
      "19108      1623     17016               4\n",
      "19109      1623     14701               3\n",
      "19110      1623     17473               4\n",
      "\n",
      "[19111 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     14871               3\n",
      "1             0     13288               3\n",
      "2             0     10231               1\n",
      "3             0      9414               1\n",
      "4             0     15902               3\n",
      "...         ...       ...             ...\n",
      "42586      5925     11262               2\n",
      "42587      5925     20132               5\n",
      "42588      5925     25023               6\n",
      "42589      5925     25281               6\n",
      "42590      5925      1949               0\n",
      "\n",
      "[42591 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     11739               2\n",
      "1             0     10542               2\n",
      "2             0     10994               2\n",
      "3             0     24497               6\n",
      "4             0     11944               2\n",
      "...         ...       ...             ...\n",
      "19417      1303     17059               4\n",
      "19418      1303     10425               2\n",
      "19419      1303     17715               4\n",
      "19420      1303     24902               6\n",
      "19421      1303      2213               0\n",
      "\n",
      "[19422 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     13118               3\n",
      "1             0     13153               3\n",
      "2             0     24554               6\n",
      "3             0     15588               3\n",
      "4             0     15103               3\n",
      "...         ...       ...             ...\n",
      "52423      5902     29613               7\n",
      "52424      5902     19147               5\n",
      "52425      5902     13676               3\n",
      "52426      5902     20235               5\n",
      "52427      5902     15885               3\n",
      "\n",
      "[52428 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     13227               3\n",
      "1             0      6338               0\n",
      "2             0      8546               1\n",
      "3             0     12402               3\n",
      "4             0     14854               3\n",
      "...         ...       ...             ...\n",
      "49277      4533      9837               1\n",
      "49278      4533     20415               5\n",
      "49279      4533     30256               7\n",
      "49280      4533     38164               9\n",
      "49281      4533     12087               2\n",
      "\n",
      "[49282 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     24365               6\n",
      "1             0     25026               6\n",
      "2             0     24166               6\n",
      "3             0     39561              10\n",
      "4             0     31833               8\n",
      "...         ...       ...             ...\n",
      "15422       994     18729               5\n",
      "15423       994      8783               1\n",
      "15424       994     18542               5\n",
      "15425       994      8596               1\n",
      "15426       994     24459               6\n",
      "\n",
      "[15427 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     28738               7\n",
      "1             1     14206               3\n",
      "2             1     22047               5\n",
      "3             1     15720               3\n",
      "4             1     15264               3\n",
      "...         ...       ...             ...\n",
      "13531      1017     17445               4\n",
      "13532      1017     31492               8\n",
      "13533      1017     32332               8\n",
      "13534      1017     10962               2\n",
      "13535      1018     31975               8\n",
      "\n",
      "[13536 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     38497              10\n",
      "1             0     41042              10\n",
      "2             0     39070              10\n",
      "3             0     39600              10\n",
      "4             0     30734               7\n",
      "...         ...       ...             ...\n",
      "25628      1769     24526               6\n",
      "25629      1769     25021               6\n",
      "25630      1770     16693               4\n",
      "25631      1770     42715              10\n",
      "25632      1770     43039              11\n",
      "\n",
      "[25633 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     14757               3\n",
      "1             0      9265               1\n",
      "2             0     24867               6\n",
      "3             0     12791               3\n",
      "4             0     15260               3\n",
      "...         ...       ...             ...\n",
      "33071      2211     27141               7\n",
      "33072      2211     27991               7\n",
      "33073      2211     26987               7\n",
      "33074      2211     14221               3\n",
      "33075      2211     24647               6\n",
      "\n",
      "[33076 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     31771               8\n",
      "1             0     16410               4\n",
      "2             0     11069               2\n",
      "3             0     17288               4\n",
      "4             0     16853               4\n",
      "...         ...       ...             ...\n",
      "28474      1656     43821              11\n",
      "28475      1656     39983              10\n",
      "28476      1656     43436              11\n",
      "28477      1656      6150               0\n",
      "28478      1656     34213               9\n",
      "\n",
      "[28479 rows x 3 columns],        new cell  ref cell  original batch\n",
      "0             0     11431               2\n",
      "1             0     12005               2\n",
      "2             0     32366               8\n",
      "3             0     11050               2\n",
      "4             0     11530               2\n",
      "...         ...       ...             ...\n",
      "44156      3069     28583               7\n",
      "44157      3069     37411               9\n",
      "44158      3069     25310               7\n",
      "44159      3069      6572               0\n",
      "44160      3069     24359               6\n",
      "\n",
      "[44161 rows x 3 columns]], None)\n"
     ]
    }
   ],
   "source": [
    "cdata = sc.external.pp.mnn_correct(looms[0],\n",
    "                                   looms[1],\n",
    "                                   looms[2],\n",
    "                                   looms[3],\n",
    "                                   looms[4],\n",
    "                                   looms[5],\n",
    "                                   looms[6],\n",
    "                                   looms[7],\n",
    "                                   looms[8],\n",
    "                                   looms[9],\n",
    "                                   looms[10],\n",
    "                                   looms[11],\n",
    "                                   looms[12],\n",
    "                                   looms[13],\n",
    "                                   looms[14],\n",
    "                                   looms[15],\n",
    "                                   looms[16],\n",
    "                                   svd_dim=50, batch_key=\"Sample\", save_raw = True, var_subset = hvars)\n",
    "print(cdata)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "0f030926",
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "17\n",
      "AnnData object with n_obs × n_vars = 56769 × 28002\n",
      "    obs: 'Sample', 'Timepoint', 'Age', 'Batch', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'batch'\n",
      "    var: 'Accession', 'Chromosome', 'End', 'Start', 'Strand', 'mt', 'n_cells_by_counts-0', 'mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'highly_variable-0', 'means-0', 'dispersions-0', 'dispersions_norm-0', 'mean-0', 'std-0', 'n_cells_by_counts-1', 'mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'highly_variable-1', 'means-1', 'dispersions-1', 'dispersions_norm-1', 'mean-1', 'std-1', 'n_cells_by_counts-10', 'mean_counts-10', 'pct_dropout_by_counts-10', 'total_counts-10', 'highly_variable-10', 'means-10', 'dispersions-10', 'dispersions_norm-10', 'mean-10', 'std-10', 'n_cells_by_counts-11', 'mean_counts-11', 'pct_dropout_by_counts-11', 'total_counts-11', 'highly_variable-11', 'means-11', 'dispersions-11', 'dispersions_norm-11', 'mean-11', 'std-11', 'n_cells_by_counts-12', 'mean_counts-12', 'pct_dropout_by_counts-12', 'total_counts-12', 'highly_variable-12', 'means-12', 'dispersions-12', 'dispersions_norm-12', 'mean-12', 'std-12', 'n_cells_by_counts-13', 'mean_counts-13', 'pct_dropout_by_counts-13', 'total_counts-13', 'highly_variable-13', 'means-13', 'dispersions-13', 'dispersions_norm-13', 'mean-13', 'std-13', 'n_cells_by_counts-14', 'mean_counts-14', 'pct_dropout_by_counts-14', 'total_counts-14', 'highly_variable-14', 'means-14', 'dispersions-14', 'dispersions_norm-14', 'mean-14', 'std-14', 'n_cells_by_counts-15', 'mean_counts-15', 'pct_dropout_by_counts-15', 'total_counts-15', 'highly_variable-15', 'means-15', 'dispersions-15', 'dispersions_norm-15', 'mean-15', 'std-15', 'n_cells_by_counts-16', 'mean_counts-16', 'pct_dropout_by_counts-16', 'total_counts-16', 'highly_variable-16', 'means-16', 'dispersions-16', 'dispersions_norm-16', 'mean-16', 'std-16', 'n_cells_by_counts-17', 'mean_counts-17', 'pct_dropout_by_counts-17', 'total_counts-17', 'highly_variable-17', 'means-17', 'dispersions-17', 'dispersions_norm-17', 'mean-17', 'std-17', 'n_cells_by_counts-2', 'mean_counts-2', 'pct_dropout_by_counts-2', 'total_counts-2', 'highly_variable-2', 'means-2', 'dispersions-2', 'dispersions_norm-2', 'mean-2', 'std-2', 'n_cells_by_counts-3', 'mean_counts-3', 'pct_dropout_by_counts-3', 'total_counts-3', 'highly_variable-3', 'means-3', 'dispersions-3', 'dispersions_norm-3', 'mean-3', 'std-3', 'n_cells_by_counts-4', 'mean_counts-4', 'pct_dropout_by_counts-4', 'total_counts-4', 'highly_variable-4', 'means-4', 'dispersions-4', 'dispersions_norm-4', 'mean-4', 'std-4', 'n_cells_by_counts-5', 'mean_counts-5', 'pct_dropout_by_counts-5', 'total_counts-5', 'highly_variable-5', 'means-5', 'dispersions-5', 'dispersions_norm-5', 'mean-5', 'std-5', 'n_cells_by_counts-6', 'mean_counts-6', 'pct_dropout_by_counts-6', 'total_counts-6', 'highly_variable-6', 'means-6', 'dispersions-6', 'dispersions_norm-6', 'mean-6', 'std-6', 'n_cells_by_counts-7', 'mean_counts-7', 'pct_dropout_by_counts-7', 'total_counts-7', 'highly_variable-7', 'means-7', 'dispersions-7', 'dispersions_norm-7', 'mean-7', 'std-7', 'n_cells_by_counts-8', 'mean_counts-8', 'pct_dropout_by_counts-8', 'total_counts-8', 'highly_variable-8', 'means-8', 'dispersions-8', 'dispersions_norm-8', 'mean-8', 'std-8', 'n_cells_by_counts-9', 'mean_counts-9', 'pct_dropout_by_counts-9', 'total_counts-9', 'highly_variable-9', 'means-9', 'dispersions-9', 'dispersions_norm-9', 'mean-9', 'std-9'\n",
      "    layers: 'ambiguous', 'matrix', 'spliced', 'unspliced'\n"
     ]
    }
   ],
   "source": [
    "print(len(processed_adatas))\n",
    "adata = adata.concatenate(processed_adatas, join = \"outer\")\n",
    "print(adata)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "f68e6652",
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "(53699, 28002)"
      ]
     },
     "execution_count": 9,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "corr_data = cdata[0][:,hvars]\n",
    "corr_data.X.shape"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "id": "b5a29b0e",
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "... storing 'Timepoint' as categorical\n",
      "... storing 'Batch' as categorical\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<Figure size 1374x800 with 4 Axes>"
      ]
     },
     "metadata": {
      "image/png": {
       "height": 717,
       "width": 1194
      }
     },
     "output_type": "display_data"
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "2022-03-09 13:57:51.984341: W tensorflow/stream_executor/platform/default/dso_loader.cc:64] Could not load dynamic library 'libcudart.so.11.0'; dlerror: libcudart.so.11.0: cannot open shared object file: No such file or directory\n",
      "2022-03-09 13:57:51.984392: I tensorflow/stream_executor/cuda/cudart_stub.cc:29] Ignore above cudart dlerror if you do not have a GPU set up on your machine.\n"
     ]
    }
   ],
   "source": [
    "# the variable genes defined are used by default by the pca function, \n",
    "# now we want to run on all the genes in the dataset\n",
    "sc.tl.pca(corr_data, svd_solver = 'arpack', use_highly_variable = False)\n",
    "sc.pl.pca(corr_data, components = ['1,2','3,4','5,6','7,8'], ncols=2, color='Sample')\n",
    "\n",
    "\n",
    "# tSNE\n",
    "sc.tl.tsne(corr_data, n_pcs = 50)\n",
    "# UMAP, first with neighbor calculation \n",
    "sc.pp.neighbors(corr_data, n_pcs = 50, n_neighbors = 20)\n",
    "sc.tl.umap(corr_data)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "91f9813d",
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "<AxesSubplot:title={'center':'MNN Corrected tsne'}, xlabel='tSNE1', ylabel='tSNE2'>"
      ]
     },
     "execution_count": 17,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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AAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsA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AAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKBwn8WVmXn2zLxhZj40M3fMzE0z846Z+Y6ZecgO537JzPzjmfmVmbl5Zu6cmf81M++dmVfPzNfOzEW7mg8AAAAAABwuF+96wMxckuTHkjz7hI8eufn70iQvmZmvW2u98xzOfUiSf5rkG5PMCR9/zubvjyb5S0keluTWczUbAAAAAAA4vHYaVzZPhLwhydM2pz6S5Nok1ye5NMlVSa5McjTJdTNz5Vrr187B3EuT/EySJ2xOHU/yb5N8IMnHkzw4yR9O8tTsBRYAAAAAAIB92fWTKy/MvWHl+iRfudb6yNbnPzQz35/kZdl7euQ1SZ5yDua+PveGlX+U5G+vte48yb7vnJnPSXL7OZgJAAAAAABcAHb2zpXNUysv3zr1/BPCyj2uSfL+zfrJM/NVZzn36iRfvTl89VrrmlOElSTJWut/r7XuOpuZAAAAAADAhWOXL7R/SpLLN+u3rbXed7JNa627k7xq69RVZzn3ms2/tyf5zrO8FgAAAAAAwO+xy7jy9K31dWfY++ZTfK8yM1cm+fzN4b9ba9120GsBAAAAAACczC7jyuO21u853ca11oeTHNscPmpmHnnAmV++tX5XkszMn52Z62bmwzNz58z875n56Zn5ppnZ9TtnAAAAAACAQ2aXceXRW+sb97F/e8+jT7nr9J6wtb5pZt6Y5I3ZexrmUUk+PXs/VfaMJD+S5H0zc8UBZwEAAAAAABegXT658dCt9S372P/RU3y3cfnW+hXZizR3Jnldkl9KcleSL0rywiSXZu/pml+YmcevtT7WDJqZI2fYcllzPQAAAAAA4P5hl3Hlkq31HfvY/4mt9YMPOPNhW+tHZy/Y/Km11ge2zr9+Zv5Jkp9L8gVJPjfJ9yb5tnLWsTNvAQAAAAAADptd/izY+XDi/+dlJ4SVJL/7jpev3zp19cw8ZKd3BgAAAAAAHAq7fHLl9tz7JMkDN8en86Ct9W0HnLn9vduT/PipNq61PjAz70zyx5N8RpIrk7y5mHX0DJ9fluQ9xfUAAAAAAID7gV3GlVtzb1x5RM4cVx5+wncP4je31h9ca33yDPvfm724kiSf1wxaax0/3ecz01wOAAAAAAC4n9jlz4LdsLW+Yh/7t/fccMpdp/dft9Yf38f+7T1+FgwAAAAAADijXcaVD26tn3i6jTPzqNz7M1s3rbVuPuDM7fer/L597N/es58YAwAAAAAAXOB2GVfesrV++hn2PmNrfd1ZzHxzkrVZP25mHnCG/U/YWh/0aRkAAAAAAOACssu48rYkH96sv2JmHn+yTTNzUZKXbp36iYMO3LwH5W2bw0uSfP2p9s7MF+Xe963cluTtB50LAAAAAABcOHYWV9Zadyd5xdap183MZ59k6/cl+eLN+u1rrZ852fVm5uqZWZu/t55m9HdurX9gZr7wJNd6VJIf2zr1qrXWJ05zTQAAAAAAgCTJxTu+/rVJnpvkqUkek+QDM3NtkuuTXJrkqiRfttl7a5IXne3AtdZ/npl/mOSaJA9P8u6Z+dEkv5TkruyFnBdu5ifJe5N8z9nOBQAAAAAALgw7jStrrbtm5muSvD7JM5NcluS7T7L1eJLnrbV+9RzN/Rszc3f2AstnJPnWzd+JfjbJn1tr3XEu5gIAAAAAAIffLt+5kiRZa9221npWkuckeVOSY0nuTHJLkndlL4A8dq31jnM8928l+aNJfjDJf83ee1XuSPI/s/delz+z1vrqtdZvnsu5AAAAAADA4TZrrfN9D4fSzBzJXkjKsWPHcuTIkfN8RwAAAAAA8Knl+PHjOXr06D2HR9dax8/n/ezXzp9cAQAAAAAAOEzEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIK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AAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAoiCsAAAAAAAAFcQUAAAAAAKAgrgAAAAAAABTEFQAAAAAAgIK4AgAAAAAAUBBXAAAAAAAACuIKAAAAAABAQVwBAAAAAAAo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uAIAAAAAANAgrgAAAAAAADSIKwAAAAAAAA3iCgAAAAAAQIO4AgAAAAAA0CCuAAAAAAAANIgrAAAAAAAADeIKAAAAAABAg7gCAAAAAADQIK4AAAAAAAA0iCsAAAAAAAAN4goAAAAAAECDuAIAAAAAANAgrgAAAAAAADSIKwAAAAAAAA3iCgAAAAAAQIO4AgAAAAAA0CCuAAAAAAAANIgrAAAAAAAADeIKAAAAAABAg7gCAAAAAADQIK4AAAAAAAA0iCsAAAAAAAAN4goAAAAAAECDuAIAAAAAANAgrgAAAAAAADSIKwAAAAAAAA3iCgAAAAAAQIO4AgAAAAAA0CCuAAAAAAAANIgrAAAAAAAADeIKAAAAAABAg7gCAAAAAADQIK4AAAAAAAA0iCsAAAAAAAAN4goAAAAAAECDuAIAAAAAANAgrgAAAAAAADSIKwAAAAAAAA3iCgAAAAAAQIO4AgAAAAAA0CCuAAAAAAAANIgrAAAAAAAADeIKAAAAAABAg7gCAAAAAADQIK4AAAAAAAA0iCsAAAAAAAAN4goAAAAAAECDuAIAAAAAANAgrgAAAAAAADSIKwAAAAAAAA3iCgAAAAAAQIO4AgAAAAAA0CCuAAAAAAAANIgrAAAAAAAADeIKAAAAAABAg7gCAAAAAADQIK4AAAAAAAA0iCsAAAAAAAAN4goAAAAAAECDuAIAAAAAANAgrgAAAAAAADSIKwAAAAAAAA3iCgAAAAAAQIO4AgAAAAAA0CCuAAAAAAAANIgrAAAAAAAADeIKAAAAAABAg7gCAAAAAADQIK4AAAAAAAA0iCsAAAAAAAAN4goAAAAAAECDuAIAAAAAANAgrgAAAAAAADSIKwAAAAAAAA3iCgAAAAAAQIO4AgAAAAAA0CCuAAAAAAAANIgrAAAAAAAADeIKAAAAAABAg7gCAAAAAADQIK4AAAAAAAA0iCsAAAAAAAAN4goAAAAAAECDuAIAAAAAANAgrgAAAAAAADSIKwAAAAAAAA3iCgAAAAAAQIO4AgAAAAAA0CCuAAAAAAAANIgrAAAAAAAADeIKAAAAAABAg7gCAAAAAADQIK4AAAAAAAA0iCsAAAAAAAAN4goAAAAAAECDuAIAAAAAANAgrgAAAAAAADSIKwAAAAAAAA3iCgAAAAAAQIO4AgAAAAAA0CCuAAAAAAAANIgrAAAAAAAA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\n",
      "text/plain": [
       "<Figure size 800x800 with 1 Axes>"
      ]
     },
     "metadata": {
      "image/png": {
       "height": 811,
       "width": 811
      },
      "needs_background": "light"
     },
     "output_type": "display_data"
    },
    {
     "data": {
      "text/plain": [
       "<Figure size 600x400 with 1 Axes>"
      ]
     },
     "metadata": {
      "image/png": {
       "height": 365,
       "width": 615
      }
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "import matplotlib.pyplot as plt\n",
    "\n",
    "fig, axs = plt.subplots(1, 1, figsize=(8,8),constrained_layout=True)\n",
    "sc.pl.tsne(corr_data, color=\"Sample\", title=\"MNN Corrected tsne\", show=False)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "id": "118288d5",
   "metadata": {},
   "outputs": [],
   "source": [
    "save_file = './scanpy_mnn_corrected.h5ad'\n",
    "corr_data.write_h5ad(save_file)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.10.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}