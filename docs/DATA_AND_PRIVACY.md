# Data and Privacy

## What is committed

- Small model checkpoints and prediction tables.
- Public-data experiment summaries.
- Submission and provenance manifests.
- Report source, figures, PDFs, and methods forms.

## What is excluded

- Patient FASTQ files.
- The complete supplied child VCF and phenotype document.
- Local OpenCRAVAT databases.
- Hugging Face tokens or other credentials.
- Multi-gigabyte public benchmark downloads.

The gated challenge bundle is approximately 85 GB compressed and must remain on approved controlled storage. The pipeline performs annotation locally and does not send patient variants to public annotation APIs.

## Local layout

Set `RDKID_DATA_DIR` to keep large inputs outside the checkout:

```bash
export RDKID_DATA_DIR=/secure/path/rare_disease_data
rdkid doctor --profile wgs
```

The repository-local `data/open/` and `data/private/` locations are ignored as a fallback. Public datasets should retain their upstream filenames so their hashes can be checked against [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

## Safe logging

Do not place tokens in command history, manifests, notebooks, or issue reports. The downloader reads supported token environment variables or the Hugging Face CLI credential store and masks tokens in diagnostic output.

