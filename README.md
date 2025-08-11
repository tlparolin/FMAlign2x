# FMAlign2x - An extended version of FMAlign2 for aligning multiple ultra-long sequences

FMAlign2x is an extended version of [FMAlign2](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btae014/7515251?searchresult=1) that enables in-memory alignment of segments found between MEMs using the SPOA library (SIMD partial order alignment tool). This feature aims to reduce the computational load of the selected primary alignment method (MAFFT, HAlign2, or HAlign3).

## Table of Contents

- [Installation](#Installation)
- [Usage](#Usage)
- [Data](#Data)
- [Issue](#Issue)
- [Uninstall](#Uninstall)
- [Related](#Related)
- [Citation](#Citation)
- [License](#License)

## Installation

We recommend using [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/) to manage the build environment for FMAlign2x.

This ensures that all dependencies are correctly isolated and reproducible.

### Requirements

- **Linux** (See note below)
- **CMake**
- **GNU Make**
- **g++** (version with OpenMP and C++23 support)

> Note: FMAlign2x can be built inside WSL, but native Linux is recommended for better performance and reliability. Compilation on WSL has not been extensively tested.

---

### Conda-based Installation

1. **Clone the repository**

   ```
   git clone https://github.com/tlparolin/FMAlign2x.git
   cd FMAlign2x
   ```

2. **Create the Conda environment**

   ```
   conda env create --file environment.yml
   conda activate fmalign2x
   ```

3. **Build the project**

   ```
   mkdir build
   cd build
   cmake -DCMAKE_BUILD_TYPE=Release [-DM64=ON] ..
   cmake --build . -j
   ```

   > Note: We provide two compilation modes: 32-bit and 64-bit. In most cases, the 32-bit mode is sufficient to handle most data. However, if the concatenated length of all sequences exceeds the range of uint32_t (4294967295), you should add the -DM64=ON parameter when compiling the program to generate a 64-bit executable.

4. **Run the program**

   ./FMAlign2 [options]

5. **Optional Features**

   Use -x 0 to disable in-memory alignment of small blocks between MEMs.

---

**Please note that if you choose halign2 and halign3 as your multiple sequence alignment methods, make sure you have Java environment installed.** To check the version of Java installed on your system, you can open a command prompt or terminal and execute the following command:

```bash
java -version
```

This will display the installed Java version information.

If you don't have Java installed or if the installed version is not compatible, you can follow these steps to install Java on Linux:

1. Update Package Lists: Run the command `sudo apt update` to update the package lists on your system.
2. Install OpenJDK: Run the command `sudo apt install default-jdk` to install the default version of OpenJDK.
3. Verify Installation: After the installation is complete, run `java -version` to verify that Java is installed and the correct version is displayed.

Once you have Java installed and verified the version, you should be able to use halign2 and halign3 for multiple sequence alignment.

## Usage

> Reminder: Please ensure that all external files (such as MAFFT, HALIGN, etc.) are properly copied to their corresponding directories. Pay close attention to the relative paths between FMAlign2 and the ext folder to avoid issues during execution.

```shell
./FMAlign2x -i /path/to/data [other options]
```

if you want to show the parameters details:

```
./FMAlign2 -h
```

Parameters Details:

- -i [file path] **[required]** The path to the input file.
- -o [output path] [default: ouput.fmaligned2.fasta] The path to the output file.
- -p [package] [default: mafft] The MSA method used in parallel align. for example, [**halign3**](https://github.com/malabz/HAlign-3), [**halign2**](https://github.com/ShixiangWan/HAlign2.0) and [**mafft**](https://mafft.cbrc.jp/alignment/software/).
- -t [int] [default: cpu number] The maximum number of threads that the program runs, the recommended setting is the number of CPUs.
- -l [int] [default: square root of mean length] The minimum length of MEMs, the default value is square root of mean length.
- -c [float] [default: 1] A floating-point parameter that specifies the minimum coverage across all sequences, with values ranging from 0 to 1.
- -f [mode] [default: global or local] The filter MEMs mode. The default setting is that if sequence number less 100, **accurate** mode otherwise **global** mode.
- -x [int] [default:1] In-memory alignment of small blocks between MEMs with SPOA.
- -d [int] [default:0] Depth of recursion, you could ignore it.
- -v [int] [default:1] Verbose option, 0 or 1. You could ignore it.
- -h [help] print help information

---

We will demonstrate with the example data `mt1x.fasta`.

```shell
./FMAlign2 -i ./data/mt1x.fasta -l 20 -c 1 -p mafft -f gloabl -o output.fmaligned2.fasta
```

This command specifies the following options:

- Input data: `mt1x.fasta` located in the `data` folder.
- Minimum length of MEMs: 20.
- Sequence coverage of MEMs: 1.
- Parallel alignment method: mafft
- Alignment mode: global mode.
- Output file: `output.fmaligned2.fasta` will be generated in the FMAlign2 directory.

After running this command, you will obtain the aligned output in the `output.fmaligned2.fasta` file.

---

If you want to evaluate the generated alignment results, you can run the `sp.py` script (requires a Python environment) with the following parameters:

```shell
python sp.py --input output.fmalign2.fasta --match 0 --mismatch 1 --gap1 2 --gap2 2
```

This command will calculate and print the SP (Sum-of-Pairs) score for the multiple sequence alignment results. The `--input` parameter specifies the input alignment file (`output.fmalign2.fasta` in this case), and the `--match`, `--mismatch`, `--gap1`, and `--gap2` parameters define the scoring scheme for matches, mismatches, and gap penalties.

By running this command, you will obtain the SP score, which provides an evaluation of the alignment quality.

## Data

Data can be assessed in [data](https://github.com/tlparolin/FMAlign2x/tree/master/data) fold. All the data is compressed using xz compression. Before using it, please decompress the files.

Here are the methods to decompress the files on different operating systems:

**Decompressing on Linux:**

1. Open the terminal.

2. Navigate to the directory where the compressed file is located.

3. Run the following command to decompress the file:

   ```
   xz -d filename.xz
   ```

   Replace `filename.xz` with the name of the file you want to decompress.

**Decompressing on Windows:**

1. Download and install an xz compression tool for Windows, such as [7-Zip](https://www.7-zip.org/) or [WinRAR](https://www.win-rar.com/).
2. Right-click on the compressed file.
3. Select "Extract to" or a similar option to decompress the file.

Please note that the decompressed files will occupy more disk space. Make sure you have enough disk space to store the uncompressed files.

If you need more data, you can visit http://lab.malab.cn/~cjt/MSA/datasets.html for more datasets.

## Issue

FMAlign2x is supported by [Scientific Computing Group Lab](https://www.ibilce.unesp.br/gcc). If you have any suggestions or feedback, we encourage you to provide them through the issue page on the project's repository. You can also reach out via email to thiago.parolin@unesp.br.

We value your input and appreciate your contribution to improving the project. Thank you for taking the time to provide feedback, and we will address your concerns as soon as possible.

## Uninstall

To completely uninstall **FMAlign2x** and clean up its build files, binaries, and the Conda environment, follow the steps below:

1. **Delete the build directory (CMake cache, object files, etc.)**

```bash
rm -rf build
```

This will remove all temporary CMake files and compiled objects.

2. **Remove the compiled binary**
   If you installed or copied the binary manually somewhere (e.g., bin/FMAlign2x or system-wide), remove it with:

```
rm -f FMAlign2x
```

Or if placed in a custom path:

```
rm -f /path/to/your/installation/FMAlign2x
```

3. **Remove the Conda environment**
   If you created a Conda environment specifically for FMAlign2x (e.g., named fmalign2x), deactivate and remove it:

```
conda deactivate
conda env remove -n fmalign2x
```

## Related

- [FMAlign2](https://github.com/metaphysicser/FMAlign2/): A novel fast multiple nucleotide sequence alignment method for ultra-long datasets
- [HAlign3](https://github.com/malabz/HAlign-3) and [HAlign2](https://github.com/ShixiangWan/HAlign2.0)
- [SPOA](https://github.com/rvaser/spoa): SIMD partial order alignment tool/library.
- [LibSAIS](https://github.com/IlyaGrebnov/libsais): linear-time construction of suffix array (SA), generalized suffix array (GSA), longest common prefix (LCP) array, permuted LCP (PLCP) array, Burrows-Wheeler transform (BWT) and inverse BWT.
- [WMSA](https://github.com/malabz/WMSA) and [WMSA2](https://github.com/malabz/WMSA2)
- [TPRA](https://github.com/malabz/TPRA): A refinement tool for ensembling different multiple sequence alignment results

## Citation

In progress.

## License

[Apache 2.0](https://github.com/tlparolin/FMAlign2x/blob/master/LICENSE).
FMAlign2x is based on the original FMAlign2 version developed by [Pinglu Zhang](https://github.com/metaphysicser).
