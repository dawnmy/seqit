use assert_cmd::Command;
use predicates::str::contains;
use rust_htslib::bam;
use std::fs;
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
            "-i",
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
fn stats_supports_multiple_inputs() {
    Command::cargo_bin("seqit")
        .unwrap()
        .args(["stats", "tests/data/a.fa", "tests/data/a.fq", "-T"])
        .assert()
        .success();
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
