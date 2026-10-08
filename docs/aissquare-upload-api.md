# AIS Square upload API record

The upload was completed through the authenticated AIS Square web API because
`aissq-explorer` does not implement publishing.

## Authentication

The web client stores an AIS Square bearer token as the `aisToken` item in the
AIS Square origin's browser LocalStorage. For this local workflow, the token is
kept in the ignored file at the repository root:

```text
config.local.yml
```

That file must never be committed or copied into a public repository.

## Upload sequence

1. `GET https://backend.aissquare.com/upload/prefix?uploadType=models`
2. `POST /upload/generateuploadurl` with the model filename, type `models`, and the returned prefix.
3. `PUT` the binary model to the returned signed object-storage URL.
4. `POST /upload/create` with model metadata, the prefix, license, README, and file-delete list.

The API calls use `Authorization: Bearer <aisToken>`. The public artifact
created for this project is:

```text
name: ICA-series-DPA3-frozen
AIS Square ID: 444
URL: https://store.aissquare.com/models/8070570f-2710-4a54-95ec-16736ca610a8/dpa3.pth
SHA-256: 39755352b3ddbeaccf6d0e2b95aa9030f7616b2c11d0fda10e31c213068fea5e
```

The upload metadata uses a generic `ICA Research` author entry and does not
publish the account email, user ID, private paths, or token. The China object
storage endpoint reset the TLS connection during upload, so the signed US
object-storage URL returned by AIS Square was used successfully.
