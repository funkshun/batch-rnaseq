# batch-rnaseq
Tool for batch identification of RNA BWT Datasets

Dependencies can be installed with `pip install -r requirements.txt`

## Usage

`./rnaseq [options] DIRECTORY`

The target directory can either be a single dataset containing a `comp_msbwt.npy` file or a directory containing several such directories.
All reports are dumped to stdout, so redirection is encouraged.

### Options
```shell
Positional Arguments
-p, --probes         The path to the file where the Variant Probes are saved
-r, --report         The type of report generated {txt, html, json} #HTML is not currently supported
-sdp, --sdp_lookup   File containing mapping between probe sdp vector and dataset names
-t, --threshold      Minimum necessary occurrences to count probe vote
```
### Notes
By default, rnaseq searches for `VariantProbes.csv` and `SDPpositions.csv` in the current directory

## JSON

Individual reports consist of the following structure:

```json
{
  "dataset": "dataset_name",
  "timing": {
    "load": 1.1,
    "vote": 1.1
  },
  "analysis": {
    "vote_count": 1000,
    "results": [
      {
        "cross": ["cross", "name"],
        "votes": 500
      }
    ]
  }
}
```

A combined report contains a JSON list of such entries.
