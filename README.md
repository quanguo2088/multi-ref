# Fast bootstrap and reliable readout using hidden references for DNA data storage

![Alt Text](./images/image.png)

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Kit Tree Diagram](#kit-tree-diagram)
- [Example of usage](#example-of-usage)
  - [1. Fast recovery (R = 2/3)](#1-fast-recovery-r--23)
  - [2. Bootstrap recovery with Type-I Reads only (R = 5/6)](#2-bootstrap-recovery-with-type-i-reads-only-r--56)
  - [3. Bootstrap recovery with Type-I + II Reads (R = 5/6)](#3-bootstrap-recovery-with-type-i--ii-reads-r--56)
  - [4. Bootstrap recovery with Type-I + II + III Reads (R = 5/6)](#4-bootstrap-recovery-with-type-i--ii--iii-reads-r--56)
- [Note](#note)
- [License](#license)

## Overview

Synthetic DNA is becoming a promising data storage medium for future large-scale data archiving. However, data readout from massive, unordered sequencing reads requires alignment based on overlapping regions and is complicated by diverse sequencing errors. We propose a multi-stage alignment and error correction strategy via multiple-fold hidden references, transforming the de novo readout into a resequencing-like workflow. The proposed scheme is compatible with various NGS platforms, maximizes the utilization of sequencing reads, and enables assembly-free data readout under low sequencing coverage. We provide code for readout pipelines under different error conditions, divided into two main parts：

1. **Fast recovery**: In low-error-rate scenarios, the pipeline identifies reads via sliding correlation to watermark reference. Bit-wise consensus rapidly generates soft-decision information for LDPC decoding.
2. **Bootstrap recovery**: In the presence of indels, the pipeline progressively identifies reads with distinct features using multiple-fold references. The forward–backward algorithm (FBA) generates indel-corrected probability information for reliable readout.

To facilitate evaluation, bootstrap recovery is divided into three workflows, each using progressively refined references and different combinations of read types. In each workflow, the FBA is applied to produce soft information for LDPC decoding:

- **Type-I Reads only**: Identified by aligning to the embedded watermark reference; typically free of indels or containing only end-position indels.
- **Type-I + Type-II Reads**: Adds a scaffold reference constructed from Type-I reads to recover Type-II reads, which typically contain internal indels.
- **Type-I + Type-II + Type-III Reads**: Uses a regenerative reference derived from decoding feedback to recover residual Type-III reads that fall into scaffold gaps.

The entire software is implemented in C and C++, with input and output files provided alongside the program. Executable calls are organized into modular shell scripts, enabling easy and flexible deployment across different Linux distributions.

We designed and synthesized four ~40 kb DNA sequences at different LDPC code rates: DNA-40.5Kb-DR (R = 1/4), DNA-40.32Kb-ER (R = 1/2), DNA-40.5Kb-EM (R = 2/3),and DNA-40.5Kb-MC (R = 5/6). We provide the corresponding data and recovery programs to support both fast and bootstrap recovery across all four code rates, enabling accurate readout under diverse error conditions.

## Requirements

**OS Requirements**

The program has been tested on the following operating systems:

- **Ubuntu 18.04.6 LTS**
- **Ubuntu 20.04.6 LTS**

**Software Requirements**

The following tools and dependencies are required:

- **C compiler**: Ensure `gcc` is installed.
- **C++ compiler**: Ensure `g++` is installed.
- **Edlib**: Sequence alignment library. Available at [https://github.com/Martinsos/edlib](https://github.com/Martinsos/edlib).
- **LDPC codes**: LDPC encoder/decoder by Radford M. Neal. Available at [https://github.com/radfordneal/LDPC-codes](https://github.com/radfordneal/LDPC-codes).

## Kit Tree Diagram

```
├── fast_recovery_BW/                             # Fast recovery
│   ├── R0.25/
│   ├── R0.5/
│   ├── R0.67/
│   └── R0.83/
│       ├── src/
│       ├── bin/
│       ├── configure/
│       ├── sequencing_data/
│       ├── build.sh
│       └── recover.sh

├── bootstrap_recovery_TypeIReads_FBA/            # Bootstrap recovery with Type-I Reads only
│   ├── R0.25/
│   ├── R0.5/
│   ├── R0.67/
│   └── R0.83/
│       ├── src/
│       ├── bin/
│       ├── configure/
│       ├── sequencing_data/
│       ├── build.sh
│       └── recover.sh

├── bootstrap_recovery_TypeI+IIReads_FBA/         # Bootstrap recovery with Type-I + Type-II Reads
│   ├── R0.25/
│   ├── R0.5/
│   ├── R0.67/
│   └── R0.83/
│       ├── src/
│       ├── bin/
│       ├── configure/
│       ├── sequencing_data/
│       ├── build.sh
│       └── recover.sh

├── bootstrap_recovery_TypeI+II+IIIReads_FBA/     # Bootstrap recovery with Type-I + II + III Reads
│   ├── R0.25/
│   ├── R0.5/
│   ├── R0.67/
│   └── R0.83/
│       ├── src/
│       ├── bin/
│       ├── configure/
│       ├── sequencing_data/
│       ├── build.sh
│       └── recover.sh

```

## Example of usage

### 1. Fast recovery (R = 2/3)

**Command:**

```bash
cd ./fast_recovery_BW/R0.67/
./build.sh
./recover.sh
```

#### [Step 1] Sliding correlation

**Input files:**

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `DNA-40.5Kb-EM_SE150.fastq` – real sequencing data with a raw error rate of 0.2%

**Output files:**

- `correlation_result.txt` – read alignment information (sequence, correlation peak, position, strand)

#### [Step 2] Bit-wise majority voting

**Input files:**

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `correlation_result.txt` – read alignment information from Step 1

**Output files:**

- `soft_info.txt` – consensus soft information with two columns: probability of "1" and "0" at each bit position

#### [Step 3] Soft-decision LDPC decoding

**Input files:**

- `soft_info.txt` – consensus soft information from Step 2

**Output files:**

- `recovery_image.jpg` – reconstructed image
- `recovery_bitstream.txt` – decoded binary stream (43,200 bits)

Fast recovery workflows for other code rates (R = 1/4, 1/2, and 5/6) are provided and follow the same structure and usage as the R = 2/3 example.

---

### 2. Bootstrap recovery with Type-I Reads only (R = 5/6)

**Command:**

```bash
cd ./bootstrap_recovery_TypeIReads_FBA/R0.83/
./build.sh
./recover.sh
```

#### [Step 1] Sliding correlation

**Input files:**

- `SequenceL81000NoPeriodOnly2ndFILE` – known watermark sequence
- `DNA-40.5Kb-MC-Sim.fastq` – simulated sequencing data with a raw error rate of 1.2% (including 0.6% indels)

**Output files:**

- `correlation_result.txt` – read alignment information including read sequence, correlation peak, alignment position relative to the watermark, and strand orientation
- `Type-I_reads.txt` – high-correlation reads with no indels or indels near the ends
- `lowthres_reads.txt` – low-correlation reads that typically carried internal indel errors

#### [Step 2] Forward-backward algorithm (FBA)

**Input files:**

- `SequenceL81000NoPeriodOnly2ndFILE` – known watermark sequence
- `Type-I_reads.txt` – high-correlation reads from Step 1

**Output files:**

- `symbol_probability.txt` – indel-corrected symbol probability (from Type-I reads)

#### [Step 3] Consensus soft information generation

**Input files:**

- `SequenceL81000NoPeriodOnly2ndFILE` – known watermark sequence
- `symbol_probability.txt` – indel-corrected symbol probability from Step 2

**Output files:**

- `soft_info.txt` – consensus soft information (from Type-I reads)

#### [Step 4] Soft-decision LDPC decoding

**Input files:**

- `soft_info.txt` – consensus soft information from Step 3

**Output files:**

- `recovery_image.jpg` – reconstructed image
- `recovery_bitstream.txt` – decoded binary stream (54,000 bits)
- `decodedCodeword.txt` – decoded full LDPC codeword including systematic and parity bits (64,800 bits)

---

### 3. Bootstrap recovery with Type-I + II Reads (R = 5/6)

**Command:**

```bash
cd ./bootstrap_recovery_TypeI+IIReads_FBA/R0.83/
./build.sh
./recover.sh
```

#### Stage 1: *(see above)*

#### Stage 2:

#### [Step 1–2] Scaffold reference generation and filtering Type-II reads

**Input files:**

- `Type-I_reads.txt` – high-correlation reads from Stage 1
- `lowthres_reads.txt` – low-correlation reads remaining from Stage 1

**Output files:**

- `TypeII_reads.txt` – Type-II reads that typically carried internal indel errors
- `scaffold_unaligned_reads.txt` – residual reads not aligned to the scaffold

#### [Step 3–5] FBA, consensus, and decoding

**Input files:**

- `SequenceL81000NoPeriodOnly2ndFILE` – known watermark sequence
- `TypeII_reads.txt` – Type-II reads from Step 2

**Output files:**

- `recovery_image.jpg` – reconstructed image
- `recovery_bitstream.txt` – decoded binary stream (54,000 bits)
- `decodedCodeword.txt` – decoded full LDPC codeword including systematic and parity bits (64,800 bits)

---

### 4. Bootstrap recovery with Type-I + II + III Reads (R = 5/6)

**Command:**

```bash
cd ./bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.83/
./build.sh
./recover.sh
```

#### Stage 1 and Stage 2: *(see above)*

#### Stage 3:

#### [Step 1–2] Regenerative reference generation and filtering Type-III reads

**Input files:**

- `SequenceL81000NoPeriodOnly2ndFILE` – known watermark sequence
- `decodedCodeword.txt` – decoded codeword from Stage 2
- `scaffold_unaligned_reads.txt` – residual reads not aligned to the scaffold from Stage 2

**Output files:**

- `TypeIII_reads.txt` – reads typically aligned to gap regions, identified during the final recovery stage

#### [Step 3–5] FBA, consensus, and decoding

**Input files:**

- `SequenceL81000NoPeriodOnly2ndFILE` – known watermark sequence
- `TypeIII_reads.txt` – Type-III reads from Step 1

**Output files:**

- `recovery_image.jpg` – reconstructed image
- `recovery_bitstream.txt` – decoded binary stream (54,000 bits)
- `decodedCodeword.txt` – decoded full LDPC codeword including systematic and parity bits (64,800 bits)

In addition, each bootstrap recovery workflow logs summary statistics for every independent experiment in:

`./bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.83/results/recovery_status.txt`

This file contains seven columns:

1. Experiment ID
2. Erasure rate
3. Substitution error rate
4. Recovery stage (0: not recovered, 1: Type-I reads only, 2: Type-I + II reads, 3: Type-I + II + III reads)
5. Number of Type-I reads
6. Number of Type-II reads
7. Number of Type-III reads

The bootstrap recovery workflow for other code rates (R = 1/4, 1/2, and 2/3) is provided and follows the same structure and usage as the R = 5/6 example.

## Note

For comprehensive information on the coding scheme in this study, please refer to the following literature:

- MacKay, D.J.C., and Neal, R.M. (1996) Near Shannon limit performance of low density parity check codes. Electronics Letters, 32(18): 1645–1646.
- Davey, M.C., and MacKay, D.J.C. (2001) Reliable communication over channels with insertions, deletions, and substitutions. IEEE Transactions on Information Theory, 47(2): 687–698.
- Chen, W.G., Han, M.Z., Zhou, J.T., Ge, Q., Wang, P.P., Zhang, X.C., Zhu, S.Y., Song, L.F., and Yuan, Y.J. (2021) An artificial chromosome for data storage. National Science Review, 8: nwab028.
- Ge, Q., Qin, R., Liu, S., Guo, Q., Han, C., and Chen, W.G. (2025) Pragmatic soft-decision data readout of encoded large DNA. Briefings in Bioinformatics, 26: bbaf102.

## License

This project is licensed under the MIT License. See the [LICENSE] file for details.
