use assert_cmd::Command;
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
        .args(["stats", "tests/data/a.fa", "--tabular"])
        .assert()
        .success();
}
