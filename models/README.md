# Model artifacts

The phonon and elastic workflows use a DPA3/DeepMD frozen model. The source tree contains no model binary in the review release. The two 10,442,104-byte files observed in the source project (`dpa3.pth` and `frozen_model.pth`) were byte-identical; their MD5 was `2bbbe45eb9c84c77a4de6003dfa3a74f`.

Before public release, publish one frozen inference artifact through the approved aissq model-upload workflow and replace the placeholder below with an immutable URL and SHA-256. Do not publish the full training checkpoint or training dataset in Git.

```text
artifact_name: dpa3-frozen
artifact_url: PENDING_REVIEW
sha256: 39755352b3ddbeaccf6d0e2b95aa9030f7616b2c11d0fda10e31c213068fea5e
intended_use: phonon and elastic inference
```

The referenced `aissq-explorer` CLI currently implements public model/dataset
listing, metadata lookup, and downloads only; it contains no login or upload
command. Publishing this artifact to AIS Square therefore requires the
platform's own authenticated upload interface or API credentials, which are
not present on this machine. No model was uploaded during this step.
