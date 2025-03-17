# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.5.0] - 2025-01-23
### Added
- [GRD-829](https://jira.oicr.on.ca/browse/GRD-829) added bcbio ensemble vcfs step
- additional output metrics
- parameters for bcbio ensembl step

## [2.4.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only)

## [2.3.2]
### Changed
- Changes to preprocessing script and the workflow: Add name injection to ensure we have VEP-compliant ids in the output vcf headers
  Note that we still support only matched calls (no tumor-only or normal-only calls). Names of samples are mandatory (get from olive)
  Although there is a partial support for tumor-only in the workflow (normalName is optional) at this point only matched calls
  are fully supported.

## [2.3.1]
### Changed
- Minor change to the way of how parameters are set, this fixes a small issue that became apparent after initial testing

## [2.3.0]
### Changed
- Modified workflow. DISCVRseq dropped, using custom consensus-marking script from Miguel Vasquez (mikisvaz@gmail.com)

## [2.2.0]
### Added
- Added step which runs post-processing script, new varmerge-scripts version as a default

## [Unrealeased - 2023-04-11
### Added
- Added post-processing script, new tag for varmerge-scripts (1.8 with a bug introduced, fixed in 1.9)

## [2.1.1] - 2023-02-15
### Changed
- Updated  the default varmerge-script version which fixes an issue with non-existent fields in vcf
- Updated RTs accordingly

## [2.1.0] - 2022-11-16
### Added
- Introducing DISCVR Seq Toolkit as a replacement for GATK's combineVariant,
- adding default modules
### Changed
- updating the default gatk version to the latest

## [Unreleased] - 2021-11-10
### Changed
- [GP-2890](https://jira.oicr.on.ca/browse/GP-2890)

## [2.0.2] - 2020-06-10
### Changed
- Migrate to Vidarr

## [2.0.1] - 2020-09-16
### Changed
- Fixing issues with inputs, added sorting step in preprocessing task

## [2.0.0]   - 2020-07-10
### Changed
- Converting to WDL, initial import
