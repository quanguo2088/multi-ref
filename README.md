# Fast bootstrap and reliable readout using hidden references for DNA data storage

![Alt Text](./images/image.png)

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Kit Tree Diagram](#kit-tree-diagram)
- [Example of usage](#example-of-usage)
  - [1. R0.67_fast_recovery](#1-r067_fast_recovery)
    - [Step 1: Sliding Correlation](#step-1-sliding-correlation)
    - [Step 2: Bit‑wise Majority Voting](#step-2-bitwise-majority-voting)
    - [Step 3: LDPC Decoding](#step-3-ldpc-decoding)
  - [2. R0.83_bootstrap_recovery](#2-r083_bootstrap_recovery)
    - [Stage 1](#stage-1)
    - [Stage 2](#stage-2)
    - [Stage 3](#stage-3)
- [Note](#note)
- [License](#license)

## Overview

Synthetic DNA is becoming a promising data storage medium for future large-scale data archiving. However, data
readout from massive, unordered sequencing reads requires alignment based on overlapping regions and is complicated by diverse sequencing errors. We propose a multi-stage alignment and error correction strategy via multiple-fold hidden references, transforming the de novo readout into a resequencing-like workflow. The proposed scheme is compatible with various NGS platforms, maximizes the utilization of sequencing reads, and enables assembly-free data readout under low sequencing coverage. We provide code for readout pipelines under different error conditions, divided into two main parts：

1. **Fast recovery**: In low-error-rate scenarios, reads are identified via sliding correlation with the watermark reference. Bit-wise consensus generates soft-decision information for fast data recovery.
2. **Bootstrap recovery**: In the presence of indels, reads with distinct features are progressively identified using multiple-fold references. The forward-backward algorithm is employed to generate indel-corrected probability information for soft-decision decoding.

The entire software is implemented in C and C++, with input and output files provided alongside the program. Executable calls are organized into modular shell scripts, enabling easy and flexible deployment across different Linux distributions.

We designed and synthesized four ~40 kb large DNA fragments using low-density parity-check codes (https://github.com/radfordneal/LDPC-codes, by Radford M. Neal) and superimposed watermarks , respectively encoding Tagore's poem *Dreams* (DNA-40.5Kb-DR, R = 1/4), the image *Emblem* (DNA-40.5Kb-EM, R = 2/3), the image *Earthrise* (DNA-40.32Kb-ER, R = 1/2) and the image *Milk Coronet* (DNA-40.5Kb-MC, R = 5/6). We provide the corresponding data and decoding programs to support both fast and bootstrap recovery at all four code rates, enabling accurate decoding under varying error conditions.

## Requirements

**OS Requirements**

The program has been tested on the following operating systems:

- **Ubuntu 18.04.6 LTS**
- **Ubuntu 20.04.6 LTS**

**Software Requirements**

The following tools and dependencies are required:

- **C Compilers:** Ensure gcc is installed.
- **C++ Compilers:** Ensure g++ is installed.
- **Edlib**: The Edlib should be available for sequence alignment. You can download it from the official repository: https://github.com/Martinsos/edlib

## Kit Tree Diagram

```
.
├── Fast recovery/                                # Fast recovery modules for four code rates
│   ├── R0.25_fast_recovery/                      # Fast recovery for R = 1/4
│   ├── R0.5_fast_recovery/                       # Fast recovery for R = 1/2
│   ├── R0.67_fast_recovery/                      # Fast recovery for R = 2/3
│   └── R0.83_fast_recovery/                      # Fast recovery for R = 5/6
│       ├── bin/                                  # Compiled binaries
│       ├── src/                                  # Source code
│       ├── configure/                            # Watermark and decoding parameter files
│       ├── sequencing_data/                      # Sequencing reads (FASTQ format)
│       ├── Compile.sh                            # Compilation script
│       └── Fast_recovery.sh                      # One-stage fast recovery using bit-wise consensus

├── Bootstrap recovery/                           # Bootstrap recovery modules for four code rates
│   ├── R0.25_bootstrap_recovery/                 # Bootstrap recovery for R = 1/4
│   ├── R0.5_bootstrap_recovery/                  # Bootstrap recovery for R = 1/2
│   ├── R0.67_bootstrap_recovery/                 # Bootstrap recovery for R = 2/3
│   └── R0.83_bootstrap_recovery/                 # Bootstrap recovery for R = 5/6
│       ├── bin/                                  # Compiled binaries
│       ├── src/                                  # Source code
│       ├── configure/                            # Watermark and decoding parameter files
│       ├── sequencing_data/                      # Sequencing reads (FASTQ format)
│       ├── Compile.sh                            # Compilation script
│       ├── R0.83_bootstrap_recovery_stage1.sh    # Stage 1: watermark reference identify Type-I reads for data recovery
│       ├── R0.83_bootstrap_recovery_stage2.sh    # Stage 2: scaffold reference identify Type-II reads for data recovery
│       ├── R0.83_bootstrap_recovery_stage3.sh    # Stage 3: regenerative reference identify Type-III reads for data recovery
│       └── Bootstrap_recovery_thread.sh          # Full pipeline combining all three recovery stages
```

## Example of usage

### 1. R0.67_fast_recovery

**Command:**

1. **Compilation**

```bash
./Compile.sh
./Fast_recovery.sh
```

##### [Step 1] Sliding Correlation

**Inputs:**

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `DNA-40.5Kb-EM_SE150.fastq` – Sequencing data

**Outputs:**

- `correlation_result.txt` – read alignment information

##### [Step 2] Bit‑wise Majority Voting

**Inputs:**

- `correlation_result.txt` – from Step 1
- `SequenceLengthALL_FILE001R0667` – known watermark sequence

**Outputs:**

- `soft_info.txt` – consensus soft information

##### [Step 3] LDPC Decoding

**Inputs:**

- `soft_info.txt` – from Step 2

**Outputs:**

- `recovery_image.jpg` – reconstructed image
- `recovery_bitstream.txt` – decoded binary bitstream with a length of 43,200 bits

Fast recovery workflows for other code rates (R = 1/4, 1/2, and 5/6) are provided and follow the same structure and usage as the R = 2/3 example.

---

### 2. R0.83_bootstrap_recovery

**Command:**

```bash
./Compile.sh
./Bootstrap_recovery_thread.sh
```

---

##### Stage 1

```bash
./R0.83_bootstrap_recovery_stage1.sh
```

##### [Step 1] Sliding Correlation

**Inputs:**

- `SequenceL81000NoPeriodOnly2ndFILE` – known watermark sequence
- `DNA-40.5Kb-MC-Sim.fastq` – simulated sequencing data

**Outputs:**

- `correlation_result.txt` – read alignment info
- `Type-I_reads.txt` – low-error reads

##### [Step 2] Forward-Backward Algorithm

**Inputs:**

- `Type-I_reads.txt` – from Step 1
- `SequenceL81000NoPeriodOnly2ndFILE` – watermark

**Outputs:**

- `symbol_probability.txt` – indel-corrected symbol probability (Stage 1)

##### [Step 3] Consensus Soft Information Generation

**Inputs:**

- `symbol_probability.txt` – from Step 2
- `SequenceL81000NoPeriodOnly2ndFILE` – watermark

**Outputs:**

- `soft_info.txt` – consensus soft info (Stage 1)

##### [Step 4] LDPC Decoding

**Inputs:**

- `soft_info.txt` – from Step 3

**Outputs:**

- `recovery_image.jpg` – recovery image
- `correctedBitStream.txt` – decoded bitstream
- `decodedCodeword.txt` – decoded codeword

---

##### Stage 2

```bash
./R0.83_bootstrap_recovery_stage2.sh
```

##### [Step 1] Majority Voting – Generate Scaffold Reference

**Inputs:**

- `decodedCodeword.txt` – from Stage 1

**Outputs:**

- `scaffold_ref.txt` – scaffold reference

##### [Step 2] Edlib Alignment for Correlation-failed Reads

**Inputs:**

- `scaffold_ref.txt` – from Step 1
- `lowthres_reads.txt` – from Stage 1

**Outputs:**

- `TypeII_reads.txt` – type-II reads

##### [Step 3] Forward-Backward Algorithm

**Inputs:**

- `TypeII_reads.txt` – from Step 2
- `SequenceL81000NoPeriodOnly2ndFILE` – watermark

**Outputs:**

- `symbol_probability.txt` – indel-corrected symbol probability (Stage 2)

##### [Step 4] Consensus Soft Information Generation

**Inputs:**

- `symbol_probability.txt` – from Step 3
- `SequenceL81000NoPeriodOnly2ndFILE` – watermark

**Outputs:**

- `soft_info.txt` – consensus soft info (Stage 2)

##### [Step 5] LDPC Decoding

**Inputs:**

- `soft_info.txt` – from Step 4

**Outputs:**

- `recovery_image.jpg` – recovery image
- `correctedBitStream.txt` – decoded bitstream
- `decodedCodeword.txt` – decoded codeword

---

##### Stage 3

```bash
./R0.83_bootstrap_recovery_stage3.sh
```

##### [Step 1] Decode Feedback to Generate Regenerative Reference

**Inputs:**

- `decodedCodeword.txt` – from Stage 2
- `remaining_reads.txt` – from Stage 2

**Outputs:**

- `TypeIII_reads.txt` – type-III reads

##### [Step 2] Forward-Backward Algorithm

**Inputs:**

- `TypeIII_reads.txt` – from Step 1
- `SequenceL81000NoPeriodOnly2ndFILE` – watermark

**Outputs:**

- `symbol_probability.txt` – indel-corrected symbol probability (Stage 3)

##### [Step 3] Consensus Soft Information Generation

**Inputs:**

- `symbol_probability.txt` – from Step 2
- `SequenceL81000NoPeriodOnly2ndFILE` – watermark

**Outputs:**

- `soft_info.txt` – consensus soft info (Stage 3)

##### [Step 4] LDPC Decoding

**Inputs:**

- `soft_info.txt` – from Step 3

**Outputs:**

- `recovery_image.jpg` – recovery image
- `correctedBitStream.txt` – decoded bitstream
- `decodedCodeword.txt` – decoded codeword

The bootstrap recovery workflow for R = 1/2 is provided and follows the same structure and usage as the R = 5/6 example.

## Note

For comprehensive information on the structural design, coding scheme, and sequencing data information of the yeast artificial chromosome in this study, please refer to the literature:

Chen, W.G., Han, M.Z., Zhou, J.T., Ge, Q., Wang, P.P., Zhang, X.C., Zhu, S.Y., Song, L.F., and Yuan, Y.J. (2021) An artificial chromosome for data storage. Natl. Sci. Rev., 8, nwab028.

## License
