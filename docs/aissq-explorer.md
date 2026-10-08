# aissq-explorer capability record

The reference implementation was inspected from the local `aissq_test`
checkout supplied during preparation. Its machine-specific parent path is not
included in this public documentation.

```text
<local-aissq-test-checkout>
```

It is an unofficial Python client for public AIS Square resources. Its README,
CLI parser, and client implementation agree on the following commands:

```text
python -m aissq list models|datasets
python -m aissq search <keyword> --type models|datasets
python -m aissq info <name> --type models|datasets
python -m aissq download <name> --type models|datasets --output <directory>
```

The client exposes public listing, search, detail, and download methods. Its
`aissq/config.py` states that public resources require no authentication. The
CLI has no `login`, `upload`, `generateuploadurl`, or `upload/create` command;
the client has no upload method. It therefore cannot publish a new model.

The client is still useful for checking a published artifact:

```bash
python -m aissq info ICA-series-DPA3-frozen --type models
python -m aissq download ICA-series-DPA3-frozen --type models --output ./downloads
```

Do not add credentials to this client package. Public browsing and downloads
are separate from authenticated AIS Square publishing.
