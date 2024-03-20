use clap::Arg;
use rayon::prelude::*;
use std::fs::File;
use std::io;
use std::io::{BufRead, Read};
use std::io::BufReader;
use std::io::{Write};
use std::iter::repeat_with;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::io::{Seek, SeekFrom};
use rayon::ThreadPoolBuilder;


fn main() {
    let matches = clap::App::new("gaf2pack")
        .version("0.1")
        .author("Sebastian V")
        .about("Converts GFA to pack format")
        .arg(
            clap::Arg::new("gfa")
                .short('g')
                .long("gfa")
                .takes_value(true)
                .about("Input GFA file")
                .required(true),
        )
        .arg(
            clap::Arg::new("alignment")
                .short('a')
                .long("alignment")
                .takes_value(true)
                .about("Input GAF file")
                .required(true),
        )
        .arg(
            clap::Arg::new("output")
                .short('o')
                .long("output")
                .takes_value(true)
                .about("Output file")
                .required(true),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .takes_value(true)
                .required(false)
                .about("Number of threads")
                .default_value("1"),
        )
        .arg(Arg::new("method")
            .short('m')
            .long("method")
            .required(false)
            .about("New method")
        )
        .get_matches();

    let threads = matches
        .value_of("threads")
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let graph_file = matches.value_of("gfa").unwrap();
    let alignment_file = matches.value_of("alignment").unwrap();
    let output_file = matches.value_of("output").unwrap();
    let new_method = matches.is_present("method");
    eprintln!("Starting");
    eprintln!("Threads: {}", threads);
    eprintln!("Graph: {}", graph_file);
    eprintln!("Alignment: {}", alignment_file);
    eprintln!("Output: {}", output_file);
    if new_method {
        split_method(graph_file, alignment_file, output_file, threads);
    } else {
        rayon_run(graph_file, alignment_file, output_file, threads);
    }
}

pub fn rayon_run(graph_file: &str, alignment_file: &str, output_file: &str, threads: usize) {
    eprintln!("Running rayon method");
    eprintln!("Parsing GFA file");
    let (node_size, index, res) = parse_gfa_file(graph_file);
    let result = Arc::new(
        repeat_with(|| AtomicUsize::new(0))
            .take(res)
            .collect::<Vec<_>>(),
    );
    eprintln!("Running packer");
    pack_wrapper(alignment_file, &node_size, &index, &result, threads).unwrap();
    eprintln!("Writing output file");
    write_pack(output_file, &node_size, &result).unwrap();
}

/// Fast parsing of a gfa file v1
fn parse_gfa_file(filename: &str) -> (Vec<usize>, Vec<usize>, usize) {
    let mut nodes: Vec<[usize; 2]> = Vec::new();
    let file = File::open(filename).expect("Failed to open file");
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('S') {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() >= 3 {
                let id = fields[1].to_string();
                let length: usize = fields[2].len();
                nodes.push([id.parse().unwrap(), length]);
            }
        }
    }
    // Sort nodes in case they are not
    nodes.sort();

    // Check if compact and incrementing
    if nodes[0][0] == 1 && nodes[nodes.len() - 1][0] == nodes.len() {
        eprintln!("Nodes are in order");
    } else {
        eprintln!("Nodes are not in order");
    }
    eprintln!("Number of nodes: {}", nodes.len());

    let (index, res) = find_change_points(&nodes);

    // Create index
    let nodes2: Vec<usize> = nodes.into_iter().map(|x| x[1]).collect();
    eprintln!("Index size: {:?}", index.len());
    eprintln!("Total size: {:?}", res);
    (nodes2, index, res)
}

/// Index
///
/// output1. Node -> Vector pointer to the start of the coverage vector
/// output2. Total size of the coverage vector
fn find_change_points(sorted_vec: &[[usize; 2]]) -> (Vec<usize>, usize) {
    let mut change_points = Vec::new();
    let mut c = 0;
    // If the vector is empty, return an empty list

    for x in sorted_vec.iter() {
        change_points.push(c);
        for _y in 0..x[1] {
            c += 1
        }
    }
    // Return the list of change points
    (change_points, c)
}

/// Cigar operation
#[derive(Debug, PartialEq)]
enum CigarOperation {
    Match,
    Insertion,
    Deletion,
    SkippedRegion,
    SoftClip,
    HardClip,
    Padding,
    Mismatch,
}

/// Total cigar unit
/// - length
/// - operation
#[derive(Debug)]
struct CigarUnit {
    length: usize,
    operation: CigarOperation,
}

/// Parse the cigar string
fn parse_cigar_string(cigar_string: &str) -> Vec<CigarUnit> {
    let mut units = Vec::new();
    let mut current_length = String::new();

    for c in cigar_string.chars() {
        if c.is_ascii_digit() {
            current_length.push(c);
        } else {
            let length = current_length.parse::<usize>().unwrap();
            match c {
                'M' => units.push(CigarUnit {
                    length,
                    operation: CigarOperation::Match,
                }),
                'D' => units.push(CigarUnit {
                    length,
                    operation: CigarOperation::Deletion,
                }),
                _ => {}
            };

            current_length.clear();
        }
    }
    units
}

/// Parse the path from a gaf file
fn parse_path(s: &str) -> (Vec<usize>, Vec<bool>) {
    let mut dirs = vec!['>' == s.chars().next().unwrap()];
    let mut node_id = Vec::new();

    let mut k = String::new();
    for x in s[1..].chars() {
        if x.is_ascii_digit() {
            k.push(x);
        } else {
            if x == '<' {
                dirs.push(false)
            } else {
                dirs.push(true);
            }
            node_id.push(k.parse().unwrap());
            k.clear();
        }
    }
    node_id.push(k.parse().unwrap());

    (node_id, dirs)
}

fn parse_line(line: &str, nodes: &[usize], index_index: &[usize], res: &Arc<Vec<AtomicUsize>>) {
    // Perform parsing for each line
    // Replace this with your actual parsing logic
    let mut lsplit = line.split_whitespace();

    let subpath = lsplit.nth(5).unwrap();
    if subpath == "0" {
        return;
    }

    let (node_id, dirs) = parse_path(subpath);
    //println!("{:?}", &lsplit.clone().last());
    let path_len = lsplit.next().unwrap().parse::<usize>().unwrap();
    let offset = lsplit.next().unwrap().parse::<usize>().unwrap();
    let path_name = lsplit.next().unwrap().parse::<usize>().unwrap();
    let _path_total = lsplit.next().unwrap().parse::<usize>().unwrap();

    if path_len < path_name {
        return;
    }
    let cigar_string = lsplit.last().unwrap().split(':').nth(2).unwrap();
    let cigar = parse_cigar_string(cigar_string);

    let mut index = Vec::new();
    for (dir, node) in dirs.iter().zip(node_id.iter()) {
        let size = nodes[*node - 1];

        if *dir {
            for x in 0..size {
                index.push([*node, x])
            }
        } else {
            for x in (0..size).rev() {
                index.push([*node, x])
            }
        }
    }
    let mut totalc = 0;

    for x in cigar.iter() {
        if x.operation == CigarOperation::Match {
            for _y in 0..x.length {
                let which = index[offset + totalc];
                let index = index_index[which[0] - 1];
                res[index + which[1]].fetch_add(1, Ordering::Relaxed);
                totalc += 1;
            }
        } else {
            totalc += x.length;
        }
    }
}

fn pack_wrapper(
    filename: &str,
    nodes: &[usize],
    index_index: &[usize],
    res: &Arc<Vec<AtomicUsize>>,
    threads: usize,
) -> io::Result<()> {
    // Open the file
    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();

    pool.install(|| {
        reader.lines().par_bridge().for_each(|result| {
            if let Ok(line) = result {
                parse_line(&line, nodes, index_index, res);
            }
        })
    });

    Ok(())
}

fn write_pack(outputfile: &str, nodes: &[usize], res: &Arc<Vec<AtomicUsize>>) -> io::Result<()> {
    let file = File::create(outputfile)?;
    let mut file = io::BufWriter::new(file);
    let mut c = 0;
    writeln!(file, "seq.pos\tnode.id\tnode.offset\tcoverage").expect("Not able to write");

    for (i, node_id) in nodes.iter().enumerate() {
        for x in 0..*node_id {
            writeln!(
                file,
                "{}\t{}\t{}\t{}",
                c,
                i + 1,
                x,
                res[c].load(Ordering::Relaxed)
            )
            .expect("Not able to write");
            c += 1;
        }
    }
    Ok(())
}



//---------------------------------TESTS---------------------------------//


pub fn split_method(graph_file: &str, alignment_file: &str, output_file: &str, threads: usize){
    eprintln!("Running new method");
    eprintln!("Parsing GFA file");
    let a =  count_lines_in_file(alignment_file).unwrap();
    let (node_size, index, res) = parse_gfa_file(graph_file);
    let result = Arc::new(
        repeat_with(|| AtomicUsize::new(0))
            .take(res)
            .collect::<Vec<_>>(),
    );

    let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap(); // Limit threads to 4


    println!("Number of lines in file: {}", a);
    let b = nth_newline_position(alignment_file, a/threads);
    println!("Position of 100th newline: {:?}", b);
    pool.install(|| {
        b.par_windows(2).for_each(|x| {
            tt(alignment_file, x[0], x[1], &node_size, &index, &result);
        });
    });

    write_pack(output_file, &node_size, &result).unwrap();

}

fn count_lines_in_file(file_path: &str) -> io::Result<usize> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut line_count = 0;
    for _ in reader.lines() {
        line_count += 1;
    }

    Ok(line_count)
}


fn nth_newline_position(file_path: &str, n: usize) -> Vec<u64> {
    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);

    let mut vv = Vec::new();
    let mut bytes = 0;
    let mut pp = String::new();
    vv.push(0);
    let mut i = 0;
    for line in reader.lines(){
        let line2 = line.unwrap();
        bytes += line2.len() as u64 + 1;

        if i >= n{
            if pp != line2{
                vv.push(bytes);
                i = 0;
            }
            else {
                i += 1;
            }
        } else {
            i += 1
        }
    }
    vv
}

fn tt(filename: &str, start_byte: u64, end_byte: u64, nodes: &[usize], index_index: &[usize], res: &Arc<Vec<AtomicUsize>>) {
    let file = File::open(filename).expect("Failed to open file");
    let mut reader = BufReader::new(file);

    // Move to the start position
    reader.seek(SeekFrom::Start(start_byte)).expect("dasjkdjas");

    // Calculate the number of bytes to read
    let num_bytes_to_read = (end_byte - start_byte) as usize;
    let mut buffer = vec![0; num_bytes_to_read];
    reader.read_exact(&mut buffer).expect("dasjkdhaskd");
    let mut old_line = String::new();
    for (c, line) in buffer.lines().enumerate(){
        let line2 = line.unwrap();

        if line2 != old_line{
            parse_line(&line2, nodes, index_index, res);
            old_line = line2;
        }
    }

}