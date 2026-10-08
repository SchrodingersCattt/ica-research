# Redaction report

The review bundle intentionally excludes cluster-specific information. The source project contained absolute NAS paths, a private host/port, cloud/DLC identifiers, local environment names, and private renderer paths. These were removed or replaced with environment variables and `config.local.example.yml`.

The bundle also excludes `.aissq/`, virtual environments, full training datasets, full checkpoints, raw cluster logs, and large external software installations. The smoke test scans text-based source files for the known private path and host patterns.
