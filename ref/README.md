Reference files necessary to run MutSig are too large to be hosted on GitHub, and as such,
are located in their own Google Storage Bucket: **`gs://getzlab-pcawg-mutsig2cv_nc/`**

This bucket is **requester pays**, so you will need a valid Google Cloud billing account
in order to access files within the bucket.

To retrieve all files in the bucket (~100 GB), run the following:
```
gsutil -m -u <billing project id> cp gs://getzlab-pcawg-mutsig2cv_nc/ .
```

When retrieving specific files, take care that the relative directory structure be preserved, e.g., 
bucket files in `gs://getzlab-pcawg-mutsig2cv_nc/c65e15/*` are copied to `ref/c65e15/*`.
