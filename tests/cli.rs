use assert_cmd::Command;
use predicates::prelude::PredicateBooleanExt;
use predicates::str::contains;
use rust_htslib::bam;
use std::fs;
use std::process::Command as ProcessCommand;
use tempfile::tempdir;

#[test]
fn fq2fa_works() {
    let td = tempdir().unwrap();
    let out = td.path().join("out.fa");
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "fq2fa",
            "tests/data/a.fq",
            "-o",
            out.to_str().unwrap(),
            "--format",
            "fastq",
        ])
        .assert()
        .success();
    let text = fs::read_to_string(out).unwrap();
    assert!(text.contains(">r1"));
}

#[test]
fn fq2fa_respects_explicit_format_on_stdin_duplicate_guard() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["fq2fa", "--format", "fastq"])
        .write_stdin("@r1\nACGT\n+\n!!!!\n")
        .assert()
        .success()
        .stdout(contains(">r1\nACGT"));
}

#[test]
fn fq2fa_respects_explicit_format_on_stdin() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["fq2fa", "--format", "fastq"])
        .write_stdin("@r1\nACGT\n+\n!!!!\n")
        .assert()
        .success()
        .stdout(contains(">r1\nACGT"));
}

#[test]
fn seq_auto_detects_format_from_stdin_content() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["seq"])
        .write_stdin("@r1\nACGT\n+\n!!!!\n")
        .assert()
        .success()
        .stdout(contains("@r1\nACGT"));
}

#[test]
fn paired_shuffle_is_deterministic() {
    let td = tempdir().unwrap();
    let o1 = td.path().join("o1.fq");
    let o2 = td.path().join("o2.fq");
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "shuffle",
            "--in1",
            "tests/data/a.fq",
            "--in2",
            "tests/data/b.fq",
            "-o",
            o1.to_str().unwrap(),
            "-O",
            o2.to_str().unwrap(),
            "-s",
            "7",
        ])
        .assert()
        .success();
    let t1 = fs::read_to_string(o1).unwrap();
    let t2 = fs::read_to_string(o2).unwrap();
    assert!(t1.contains("@r"));
    assert!(t2.contains("/2"));
}

#[test]
fn spike_single_runs() {
    let td = tempdir().unwrap();
    let out = td.path().join("spike.fa");
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "spike",
            "tests/data/a.fa",
            "-a",
            "tests/data/a.fa",
            "-o",
            out.to_str().unwrap(),
            "--format",
            "fasta",
            "-s",
            "2",
        ])
        .assert()
        .success();
    let text = fs::read_to_string(out).unwrap();
    assert!(text.matches('>').count() >= 6);
}

#[test]
fn sample_single_num_outputs_requested_count() {
    let td = tempdir().unwrap();
    let out = td.path().join("sample.fq");
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "sample",
            "tests/data/a.fq",
            "--format",
            "fastq",
            "-n",
            "2",
            "-o",
            out.to_str().unwrap(),
            "-s",
            "7",
        ])
        .assert()
        .success();
    let text = fs::read_to_string(out).unwrap();
    assert_eq!(text.matches('@').count(), 2);
}

#[test]
fn sample_auto_detects_format_from_stdin_content() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["sample", "-n", "1", "-s", "7"])
        .write_stdin(
            "@r1
ACGT
+
!!!!
@r2
TGCA
+
####
",
        )
        .assert()
        .success()
        .stdout(contains("@"));
}

#[test]
fn sample_paired_num_streaming_runs() {
    let td = tempdir().unwrap();
    let o1 = td.path().join("sample1.fq");
    let o2 = td.path().join("sample2.fq");
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "sample",
            "--in1",
            "tests/data/a.fq",
            "--in2",
            "tests/data/b.fq",
            "-n",
            "2",
            "-s",
            "7",
            "-o",
            o1.to_str().unwrap(),
            "-O",
            o2.to_str().unwrap(),
        ])
        .assert()
        .success();
    let t1 = fs::read_to_string(o1).unwrap();
    let t2 = fs::read_to_string(o2).unwrap();
    assert_eq!(t1.matches('@').count(), 2);
    assert_eq!(t2.matches('@').count(), 2);
}

#[test]
fn sample_rejects_partial_paired_flags() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["sample", "-i", "tests/data/a.fq", "-n", "1"])
        .assert()
        .failure()
        .stderr(contains("paired-end mode requires both read files"));
}

#[test]
fn spike_supports_paired_short_i_i() {
    let td = tempdir().unwrap();
    let o1 = td.path().join("spike1.fq");
    let o2 = td.path().join("spike2.fq");

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "spike",
            "-i",
            "tests/data/a.fq",
            "-I",
            "tests/data/b.fq",
            "-a",
            "tests/data/a.fq",
            "-A",
            "tests/data/b.fq",
            "-o",
            o1.to_str().unwrap(),
            "-O",
            o2.to_str().unwrap(),
            "-s",
            "3",
        ])
        .assert()
        .success();

    let t1 = fs::read_to_string(o1).unwrap();
    let t2 = fs::read_to_string(o2).unwrap();
    assert!(t1.contains("@r1"));
    assert!(t2.contains("/2"));
}

#[test]
fn spike_supports_single_positional_input() {
    let td = tempdir().unwrap();
    let out = td.path().join("spiked.fa");
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "spike",
            "tests/data/a.fa",
            "-a",
            "tests/data/a.fa",
            "-o",
            out.to_str().unwrap(),
            "-s",
            "5",
        ])
        .assert()
        .success();
    let t = fs::read_to_string(out).unwrap();
    assert!(t.contains(">r1"));
}

#[test]
fn spike_rejects_mixed_single_and_paired_inputs() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "spike",
            "tests/data/a.fq",
            "-i",
            "tests/data/a.fq",
            "-I",
            "tests/data/b.fq",
            "-a",
            "tests/data/a.fq",
            "-A",
            "tests/data/b.fq",
        ])
        .assert()
        .failure()
        .stderr(contains("do not mix positional single-end INPUT"));
}

#[test]
fn stats_tabular_runs() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", "tests/data/a.fa", "-T"])
        .assert()
        .success();
}

#[test]
fn stats_pretty_table_runs() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", "tests/data/a.fa"])
        .assert()
        .success();
}

#[test]
fn stats_auto_detects_format_from_stdin_content() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", "-", "-T"])
        .write_stdin("@r1\nACGT\n+\n!!!!\n")
        .assert()
        .success()
        .stdout(contains("fastq"));
}

#[test]
fn stats_supports_multiple_inputs() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", "tests/data/a.fa", "tests/data/a.fq", "-T"])
        .assert()
        .success();
}

#[test]
fn stats_reports_sequence_type_dna_rna_protein() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", "-", "-T"])
        .write_stdin(">d1\nACGTNN\n")
        .assert()
        .success()
        .stdout(contains("\tDNA\t"));

    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", "-", "-T"])
        .write_stdin(">r1\nACGUNN\n")
        .assert()
        .success()
        .stdout(contains("\tRNA\t"));

    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", "-", "-T"])
        .write_stdin(">p1\nMKWVTFISLLLLFSSAYS\n")
        .assert()
        .success()
        .stdout(contains("\tprotein\t"));
}

#[test]
fn stats_all_mode_hides_n50_l50_for_fastq_only_input() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", "tests/data/a.fq", "-T", "-a"])
        .assert()
        .success()
        .stdout(contains("q10_pct"))
        .stdout(contains("q20_pct"))
        .stdout(contains("q30_pct"))
        .stdout(predicates::str::is_match("n50\\tl50").unwrap().not());
}

#[test]
fn stats_all_mode_shows_na_for_mixed_assembly_and_read_only_metrics() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", "tests/data/a.fa", "tests/data/a.fq", "-T", "-a"])
        .assert()
        .success()
        .stdout(contains("n50\tl50"))
        .stdout(contains("q10_pct\tq20_pct\tq30_pct"))
        .stdout(contains("fastq"))
        .stdout(contains("\tNA\tNA\t"))
        .stdout(contains("fasta"));
}

#[test]
fn stats_supports_bam_input() {
    let td = tempdir().unwrap();
    let bam_path = td.path().join("tiny.bam");

    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 1000),
    );
    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    let mut rec = bam::Record::new();
    let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(4)]);
    rec.set(b"r1", Some(&cigar), b"ACGT", &[30, 30, 30, 30]);
    writer.write(&rec).unwrap();
    drop(writer);

    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", bam_path.to_str().unwrap(), "-T"])
        .assert()
        .success()
        .stdout(contains("bam"));
}

#[test]
fn seq_supports_quality_filters_and_only_id() {
    let td = tempdir().unwrap();
    let out = td.path().join("ids.txt");

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "seq",
            "tests/data/a.fq",
            "--format",
            "fastq",
            "-q",
            "40",
            "-i",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let text = fs::read_to_string(out).unwrap();
    assert_eq!(text, "r1\nr3\n");
}

#[test]
fn seq_name_defaults_to_id_and_full_name_restores_header() {
    let td = tempdir().unwrap();
    let input = td.path().join("named.fa");
    let out_ids = td.path().join("ids.txt");
    let out_full = td.path().join("full.txt");
    fs::write(&input, ">r1 alpha desc\nACGT\n>r2\nTGCA\n").unwrap();

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "seq",
            input.to_str().unwrap(),
            "--format",
            "fasta",
            "-n",
            "-o",
            out_ids.to_str().unwrap(),
        ])
        .assert()
        .success();
    assert_eq!(fs::read_to_string(&out_ids).unwrap(), "r1\nr2\n");

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "seq",
            input.to_str().unwrap(),
            "--format",
            "fasta",
            "-n",
            "--full-name",
            "-o",
            out_full.to_str().unwrap(),
        ])
        .assert()
        .success();
    assert_eq!(
        fs::read_to_string(&out_full).unwrap(),
        "r1 alpha desc\nr2\n"
    );
}

#[test]
fn seq_supports_remove_gaps_and_seq_only() {
    let td = tempdir().unwrap();
    let input = td.path().join("gapped.fa");
    let out = td.path().join("seqs.txt");
    fs::write(&input, ">r1\nAC-G T.\n").unwrap();

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "seq",
            input.to_str().unwrap(),
            "--format",
            "fasta",
            "-g",
            "-s",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let text = fs::read_to_string(out).unwrap();
    assert_eq!(text, "ACGT\n");
}

#[test]
fn head_supports_num_and_paired_i_i() {
    let td = tempdir().unwrap();
    let o1 = td.path().join("head1.fq");
    let o2 = td.path().join("head2.fq");
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "head",
            "-i",
            "tests/data/a.fq",
            "-I",
            "tests/data/b.fq",
            "-n",
            "2",
            "-o",
            o1.to_str().unwrap(),
            "-O",
            o2.to_str().unwrap(),
        ])
        .assert()
        .success();
    let t1 = fs::read_to_string(o1).unwrap();
    let t2 = fs::read_to_string(o2).unwrap();
    assert!(t1.contains("@r1"));
    assert!(!t1.contains("@r3"));
    assert!(t2.contains("@r2/2"));
}

#[test]
fn tail_supports_proportion() {
    let td = tempdir().unwrap();
    let out = td.path().join("tail.fa");
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "tail",
            "tests/data/a.fa",
            "--format",
            "fasta",
            "-p",
            "0.34",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    let text = fs::read_to_string(out).unwrap();
    assert!(text.contains(">r3"));
    assert!(!text.contains(">r1"));
}

#[test]
fn tail_proportion_streams_from_stdin() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["tail", "--format", "fastq", "-p", "0.5"])
        .write_stdin("@r1\nACGT\n+\nIIII\n@r2\nGGNN\n+\nHHHH\n@r3\nTTAA\n+\nJJJJ\n")
        .assert()
        .success()
        .stdout(contains("@r2"))
        .stdout(contains("@r3"))
        .stdout(predicates::str::contains("@r1").not());
}

#[test]
fn sort_spills_with_small_memory_limit() {
    let td = tempdir().unwrap();
    let input = td.path().join("unsorted.fa");
    let out = td.path().join("sorted.fa");
    fs::write(&input, ">r3\nTTAA\n>r1\nACGT\n>r2\nGGNN\n").unwrap();

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "sort",
            input.to_str().unwrap(),
            "--format",
            "fasta",
            "--mem",
            "1",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let text = fs::read_to_string(out).unwrap();
    assert_eq!(text, ">r1\nACGT\n>r2\nGGNN\n>r3\nTTAA\n");
}

#[test]
fn rmdup_keep_first_streams_unique_ids() {
    let td = tempdir().unwrap();
    let input = td.path().join("dup.fa");
    let out = td.path().join("dedup.fa");
    fs::write(&input, ">r1\nACGT\n>r2\nGGNN\n>r1\nTTAA\n").unwrap();

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "rmdup",
            input.to_str().unwrap(),
            "--format",
            "fasta",
            "--by",
            "id",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let text = fs::read_to_string(out).unwrap();
    assert_eq!(text, ">r1\nACGT\n>r2\nGGNN\n");
}

#[test]
fn paired_validation_reports_and_allows_override() {
    let td = tempdir().unwrap();
    let bad_r2 = td.path().join("bad2.fq");
    fs::write(
        &bad_r2,
        "@x1/2\nAAAA\n+\nIIII\n@r2/2\nCCCC\n+\nIIII\n@r3/2\nGGGG\n+\nIIII\n",
    )
    .unwrap();
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "shuffle",
            "-i",
            "tests/data/a.fq",
            "-I",
            bad_r2.to_str().unwrap(),
            "-o",
            td.path().join("o1.fq").to_str().unwrap(),
            "-O",
            td.path().join("o2.fq").to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(contains("invalid/unpaired records"));

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "shuffle",
            "-i",
            "tests/data/a.fq",
            "-I",
            bad_r2.to_str().unwrap(),
            "--allow-unpaired",
            "-o",
            td.path().join("ok1.fq").to_str().unwrap(),
            "-O",
            td.path().join("ok2.fq").to_str().unwrap(),
        ])
        .assert()
        .success();
}

#[test]
fn rename_supports_mapping_file() {
    let td = tempdir().unwrap();
    let map = td.path().join("map.tsv");
    let out = td.path().join("mapped.fa");
    fs::write(&map, "r1\talpha\nr3\tgamma\n").unwrap();

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "rename",
            "tests/data/a.fa",
            "--format",
            "fasta",
            "--map-file",
            map.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let text = fs::read_to_string(out).unwrap();
    assert!(text.contains(">alpha"));
    assert!(text.contains(">r2"));
    assert!(text.contains(">gamma"));
}

#[test]
fn rename_map_file_matches_exact_ids_only() {
    let td = tempdir().unwrap();
    let input = td.path().join("input.fa");
    let map = td.path().join("map.tsv");
    let out = td.path().join("mapped.fa");
    fs::write(
        &input,
        ">SRR23973341.1228619\nACGT\n>SRR23973341.12286190\nTGCA\n",
    )
    .unwrap();
    fs::write(&map, "SRR23973341.1228619\tTARGET\n").unwrap();

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "rename",
            input.to_str().unwrap(),
            "--format",
            "fasta",
            "--map-file",
            map.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let text = fs::read_to_string(out).unwrap();
    assert!(text.contains(">TARGET"));
    assert!(text.contains(">SRR23973341.12286190"));
    assert!(!text.contains(">TARGET0"));
}

#[test]
fn rename_supports_regex_replacement() {
    let td = tempdir().unwrap();
    let out = td.path().join("regex.fa");

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "rename",
            "tests/data/a.fa",
            "--format",
            "fasta",
            "--match-regex",
            "^r(\\d+)$",
            "--replace",
            "seq_$1",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let text = fs::read_to_string(out).unwrap();
    assert!(text.contains(">seq_1"));
    assert!(text.contains(">seq_2"));
    assert!(text.contains(">seq_3"));
}

#[test]
fn rename_help_includes_examples_and_template_placeholders() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["rename", "--help"])
        .assert()
        .success()
        .stdout(contains("Template supports prefix and index placeholders"))
        .stdout(contains("--map-file id_map.tsv"))
        .stdout(contains("--match-regex '^sample_(\\\\d+)$'"));
}

#[test]
fn seq_piped_to_head_exits_cleanly_without_broken_pipe_error() {
    let bin = assert_cmd::cargo::cargo_bin("seqit");
    let output = ProcessCommand::new("bash")
        .args([
            "-lc",
            &format!(
                "set -o pipefail; {} seq -n tests/data/a.fq | head -n 1 >/dev/null",
                bin.display()
            ),
        ])
        .output()
        .unwrap();

    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("Broken pipe"));
}

#[test]
fn grep_help_shows_numeric_paired_flags_only() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["grep", "--help"])
        .assert()
        .success()
        .stdout(contains("-1, --in1"))
        .stdout(contains("-2, --in2"))
        .stdout(contains("-i, --ignore-case"));
}

#[test]
fn rename_supports_shortcuts_for_paired_and_mode_options() {
    let td = tempdir().unwrap();
    let o1 = td.path().join("rename1.fq");
    let o2 = td.path().join("rename2.fq");
    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "rename",
            "-1",
            "tests/data/a.fq",
            "-2",
            "tests/data/b.fq",
            "-p",
            "pair_",
            "-M",
            "generate",
            "-k",
            "-o",
            o1.to_str().unwrap(),
            "-O",
            o2.to_str().unwrap(),
        ])
        .assert()
        .success();
    let t1 = fs::read_to_string(o1).unwrap();
    let t2 = fs::read_to_string(o2).unwrap();
    assert!(t1.contains("@pair_000001/1"));
    assert!(t2.contains("@pair_000001/2"));
}

#[test]
fn grep_paired_pattern_file_trims_whitespace() {
    let td = tempdir().unwrap();
    let patterns = td.path().join("ids.txt");
    let o1 = td.path().join("grep1.fq");
    let o2 = td.path().join("grep2.fq");
    fs::write(&patterns, "  r2  \n").unwrap();

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "grep",
            "-1",
            "tests/data/a.fq",
            "-2",
            "tests/data/b.fq",
            "-f",
            patterns.to_str().unwrap(),
            "-o",
            o1.to_str().unwrap(),
            "-O",
            o2.to_str().unwrap(),
        ])
        .assert()
        .success();

    let t1 = fs::read_to_string(o1).unwrap();
    let t2 = fs::read_to_string(o2).unwrap();
    assert!(t1.contains("@r2"));
    assert!(t2.contains("@r2/2"));
}

#[test]
fn grep_by_id_does_not_match_description_text() {
    let td = tempdir().unwrap();
    let input = td.path().join("desc.fastq");
    let out = td.path().join("out.fastq");
    fs::write(
        &input,
        "@x1 SRR23973341.2/1\nACGT\n+\nIIII\n@y2 nohit\nTGCA\n+\nIIII\n",
    )
    .unwrap();

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "grep",
            input.to_str().unwrap(),
            "--format",
            "fastq",
            "-p",
            "SRR23973341.2",
            "-b",
            "id",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let text = fs::read_to_string(out).unwrap();
    assert!(text.is_empty());
}

#[test]
fn grep_by_id_pattern_file_matches_exact_ids() {
    let td = tempdir().unwrap();
    let input = td.path().join("ids.fastq");
    let patterns = td.path().join("ids.txt");
    let out = td.path().join("out.fastq");
    fs::write(
        &input,
        "@SRR23973341.1228619\nACGT\n+\nIIII\n@SRR23973341.12286190\nACGT\n+\nIIII\n@SRR23973341.5050320\nACGT\n+\nIIII\n",
    )
    .unwrap();
    fs::write(&patterns, "SRR23973341.1228619\nSRR23973341.5050320\n").unwrap();

    Command::cargo_bin("seqit")
        .unwrap()
        .args([
            "grep",
            input.to_str().unwrap(),
            "--format",
            "fastq",
            "-b",
            "id",
            "-f",
            patterns.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();

    let text = fs::read_to_string(out).unwrap();
    assert!(text.contains("@SRR23973341.1228619"));
    assert!(text.contains("@SRR23973341.5050320"));
    assert!(!text.contains("@SRR23973341.12286190"));
}
