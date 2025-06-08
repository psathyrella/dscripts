### datascripts

Helper repo to simplify running partis on real data.
See the [partis manual](https://github.com/psathyrella/partis/blob/main/docs/contents.md) for all info on what that entails.

To install, clone this repo into the partis main dir, and run from there (it may require modification to work from other dirs).

General syntax is:

```./dscripts/run.py <action> --study <study> [--paired] --version <version> --base-outdir /path/to/output```

This runs the specified partis action using the yaml configuration file in `dscripts/meta/<study>/samples.yaml`.
So to run on your own data, copy the files in `dscripts/meta/test` to `dscripts/meta/<your study name>` and modify them.
For instance to run on some test files specified in [`dscripts/meta/test/samples.yaml`](https://github.com/psathyrella/dscripts/blob/main/meta/test/samples.yaml):

```./dscripts/run.py cache-parameters --study test --paired --version test-v0 --samples paired-sample-1 --base-outdir /path/to/output --dry-run --print-width 0```

Well actually with the `--dry-run` it just prints the partis commands it would run; if you're sure the commands look correct, remove `--dry-run` to actually run.
Be careful not to run this same command again, which would start a new set of jobs doing the same things as your first ones.
Once they're running, you can check status by adding `--check` (which prints the tail of the log files) or `--logfnames` (which prints the name of the log files, e.g. for piping to `|xargs less -LS`).
If you have lots of samples, you probably want to run only a few partis jobs at once, specified with `--n-max-jobs`.
By default this will start that many jobs, then wait for some to finish before starting more.
If you set `--start-n-max-and-exit`, it'll instead exit after starting the first batch of jobs.
The number of procs for each partis job is set with `--n-procs`.
To write debug output, add `--view-ascii`.

For reproducibility (and to reduce typing) you don't actually want to type these out each time; instead keep track of them in a shell script like [`dscripts/meta/test/run.sh`](https://github.com/psathyrella/dscripts/blob/main/meta/test/run.sh).

After `cache-parameters` has finished, you typically want to run partitioning, for which you'd just run the same command but with `cache-parameters` replaced with `partition`.
If you have some seed sequences (see partis manual), you can use `seed-partition` to run that action (seed sequences are specified with a `seeds.yaml` file, see example in `dscripts/meta/test/seeds.yaml`).
The `simulate` action will run partis simulation mimicking the data sample, using the parameters made by `cache-parameters`.
After running simulation, you can run other partis commands on the resulting sample by adding the dir/file as a new sample to samples.yaml with your chosen name.

There's also a group of args that run partis actions on existing output: (`--plot-partitions`, `--merge-paired-partitions`, `--get-selection-metrics`, `--infer-trees`, etc.).
Add these to a run.py partition command to run them on existing partition output that has already finished.

There are lots of other options in run.py, if they don't have a help message (or it's confusing), please open an issue to ask, this is definitely not yet a fully documented repository (sorry!).
