use rand::Rng;
use std::fs::File;

use close_file::Closable;

use serde::Deserialize;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::Write;
use std::path::PathBuf;

fn get_random() -> f32 {
    //    let t = StudentT::new(11.0).unwrap();
    //    let v = t.sample(&mut rand::thread_rng());
    let mut rng = rand::thread_rng();
    let v: f32 = rng.gen(); // generates a float between 0 and 1
    return v;
}

fn tilfaeldighed_begraenset(min: f32, max: f32) -> f32 {
    return (get_random()) * (max - min) + min;
}

pub fn make_random_genome(size: u16) -> Vec<f32> {
    let mut result: Vec<f32> = Vec::new();

    for p in 0..=(size - 1) {
        let maxval: f32 = GENE_MAX[p as usize];
        let minval: f32 = GENE_MIN[p as usize];
        let delta: f32 = maxval - minval;
        result.push((get_random() * delta) + minval);
    }
    //println!("GenomeX{:#?}", result);
    return result;
}

#[derive(Debug, PartialEq, PartialOrd, Clone)]
pub struct TsimBattery {
    pub q_cap: f32,
    pub state_of_charge: f32,
    pub v_external: f32,
    pub v_internal: f32,
    pub resistance: f32,
    pub capacitance: f32,
    pub capacitance_loss: f32,
    pub stack_current: f32,
}

impl TsimBattery {
    pub fn new(new_capacitance: f32, new_resistance: f32, new_capacitance_loss: f32) -> Self {
        TsimBattery {
            q_cap: 0.0,
            state_of_charge: 0.0,
            v_external: 0.0,
            v_internal: 0.0,
            resistance: new_resistance,
            capacitance: new_capacitance,
            capacitance_loss: new_capacitance_loss,
            stack_current: 0.0,
        }
    }
}
//#[derive(Debug, PartialEq, PartialOrd, Clone, Copy)]
type TGene = Vec<f32>;

// the model needs capacity in Farad and Resistance in Ohm
//let GENE_MAX :Vec<f32>= vec![10000.0, 100.0]; //not const since vec![] is dynamic alloc
//let GENE_MIN :Vec<f32>= vec![0.1, 0.1]; //not const since vec![] is dynamic alloc
static GENE_MAX: &'static [f32] = &[10000.0, 100.0]; //cap, res
static GENE_MIN: &'static [f32] = &[0.1, 0.1];
/// freepascal: const GENE_MIN : Array [0..1] of single = (0.1, 0.1);

const NUM_OF_GENES: u16 = 2;

#[derive(Debug, PartialEq, PartialOrd, Clone)]
pub struct TIndivid {
    pub genes: TGene,
    pub fitness: f32,
}

#[derive(Deserialize, Debug, PartialEq, Clone)]
struct DataLineRecord {
    time: f32,
    voltage: f32,
    current: f32,
}

#[derive(Debug, Default, Clone)]
pub struct TDataFrame {
    rows: Vec<DataLineRecord>,
}

pub fn read_csv(filepath: &str) -> TDataFrame {
    let file = std::fs::File::open(filepath).unwrap();
    let mut rdr = csv::Reader::from_reader(file);
    let rows: Result<Vec<_>, _> = rdr.deserialize().collect();
    let rows = rows.unwrap();
    TDataFrame { rows }
}

fn calc_voltage(
    filepath: &str,
    store_result: bool,
    sim_data: TDataFrame,
    capacitance: f32,
    resistance: f32,
) -> f32 {
    //let mut file:File ;

    let mut temp_file = PathBuf::from("unused.txt"); //Hack since append does not work the same way

    if store_result {
        temp_file = PathBuf::from(filepath); //party like it is 198x
    };
    let mut file = File::create(temp_file).unwrap();

    //let mut resultat:Vec<DataLineRecord> = vec![];
    let mut substep: u32;
    let steps_in_substeps: u32 = 4;
    let mut subdt: f32;

    let mut q_cap: f32 = 0.0;
    let mut v_cap: f32;
    let mut v_src: f32 = 0.0;
    let mut fitness: f32 = 0.0;
    let mut state_of_charge: f32;
    let mut clocktime: f32 = sim_data.rows[0].time;
    for idx in 1..sim_data.rows.len() {
        substep = 0;
        let dt: f32 = sim_data.rows[idx].time - sim_data.rows[idx - 1].time;
        subdt = dt / (steps_in_substeps as f32); // divide timestep
                                                 //        let i_charge : f32 = (sim_data.rows[idx].current + sim_data.rows[idx-1].current) * 0.5;
        let i_charge: f32 = sim_data.rows[idx].current;
        let old_i_charge: f32 = sim_data.rows[idx - 1].current;
        while substep < steps_in_substeps {
            substep = substep + 1;
            clocktime = clocktime + subdt;
            q_cap = q_cap + subdt * 0.5 * (i_charge + old_i_charge);

            if q_cap < 0.0 {
                q_cap = 0.0;
            };
            //     self.Q_cap = 0
            //     self.V_cap = 0
            //     self.I_charge = 0
            //     print("Qcap < 0")
            // else:
            v_cap = q_cap / capacitance;

            v_src = v_cap + i_charge * resistance;
            // if self.V_stack < 0:
            //     self.V_stack = 0;
            //     print("Vstack < 0 with I_charge = " + str(self.I_charge))
            //     self.I_charge = (self.V_stack-self.V_cap)/self.R_1;
            //self.P_stack_input = self.V_stack * self.I_charge;
            state_of_charge = q_cap;

            // self.Q_cap = self.Q_cap * CAP_LOSS
            // self.C = self.C * CAP_LOSS # Simulate capacity loss

            //println!("{}, {}, {}",clocktime, Q_cap, v_src);
            if store_result {
                // Write a &str in the file (ignoring the result).
                writeln!(&mut file, "{}, {}, {}", clocktime, q_cap, v_src).unwrap();
            }
        }
        let mut error: f32 = v_src - sim_data.rows[idx].voltage;
        error = error * error * dt;
        fitness = fitness + error;
    }
    //println!("Fitness = {}", fitness);
    if store_result {
        let mut _resultat = file.close(); //closing of files is handled different from most languages as a separate set of errors (see: https://docs.rs/close-file/latest/close_file/ )
    }
    return fitness;
}

impl TIndivid {
    pub fn new(genes: TGene, fitness: f32) -> Self {
        TIndivid { genes, fitness }
    }
    pub fn mutate(
        &mut self,
        current_gen: u32,
        max_gen: u32,
        mutation_intensity: f32,
        mutation_rate: f32,
    ) {
        for p in 0..self.genes.len() {
            let maxval: f32 = GENE_MAX[p as usize];
            let minval: f32 = GENE_MIN[p as usize];
            let delta: f32 = maxval - minval;
            let mut noise: f32 = 0.0;
            let tal = self.genes[p].to_owned();
            if get_random() < mutation_rate {
                noise = ((1.0 -(current_gen as f32) / (max_gen as f32))) * mutation_intensity * delta;
                    //- (delta * 0.5 as f32);
                    self.genes[p] = tilfaeldighed_begraenset(tal - noise, tal + noise);
            } else {
                self.genes[p] = tal;
            }
        }
    }
    pub fn limit_values(&mut self) {
        for p in 0..self.genes.len() {
            let maxval: f32 = GENE_MAX[p as usize];
            let minval: f32 = GENE_MIN[p as usize];
            if self.genes[p] > maxval {
                self.genes[p] = maxval;
            }
            if self.genes[p] < minval {
                self.genes[p] = minval;
            }
        }
    }
    pub fn calc_fitness(
        &mut self,
        sim_data: TDataFrame,
        filepath: &str,
        store_result: bool,
    ) -> f32 {
        let capacitance: f32 = self.genes[0];
        let resistance: f32 = self.genes[1];
        self.fitness = calc_voltage(filepath, store_result, sim_data, capacitance, resistance);
        return self.fitness;
    }
    pub fn make_random_genome(&mut self) {
        self.genes.clear();
        self.genes = make_random_genome(NUM_OF_GENES);
    }
}

/******************************************************************
    Name        : crossover
    Input       : Two TIndivid that should be crossed
    resultat    : Based on fitness the resulting genes are crossed. Fitness in TIndivid is set to -1.0 (no fitness)
    Description : Simple crossover operator for floats
******************************************************************/
fn crossover(item1: TIndivid, item2: TIndivid) -> TIndivid {
    let mut resultat = TIndivid::new(make_random_genome(NUM_OF_GENES), 25.0); //Values are not used, but we reserve space for record
                                                                              // item1.genes.to_vec(); //to_vec Copies content not just reference
    for p in 0..item1.genes.len() {
        resultat.genes[p] = (item1.genes[p] * item1.fitness + item2.genes[p] * item2.fitness)
            / (item1.fitness + item2.fitness);
    }
    resultat.fitness = -1.0; //Destroy fitness value
    return resultat;
}

#[derive(Debug, PartialEq, PartialOrd, Clone)]
pub struct TPopulation {
    pub population: Vec<TIndivid>,
    pub population_size: u32,
    pub current_generation: u32,
    pub max_generation: u32,
    pub best_fitness: Vec<f32>,
    pub avg_fitness: Vec<f32>,
    pub worst_fitness: Vec<f32>,

    /// in range 0-1
    pub mutation_rate: f32,

    ///
    pub mutation_intensity: f32,

    /// in range 0-1
    pub elite_size: f32,

    /// in range 0-1
    pub crossover_size: f32,
}

impl TPopulation {
    pub fn new(max_generation: u32, mutation_rate: f32, mutation_intensity: f32) -> Self {
        TPopulation {
            population: vec![],
            population_size: 0,
            current_generation: 0,
            max_generation: max_generation,
            best_fitness: vec![],
            avg_fitness: vec![],
            worst_fitness: vec![],
            mutation_rate: mutation_rate,
            mutation_intensity: mutation_intensity,
            elite_size: 0.1,
            crossover_size: 0.2,
        }
    }

    pub fn make_population(&mut self, population_sizer: u32) {
        for _p in 0..population_sizer {
            self.population.push(TIndivid::new(
                make_random_genome(NUM_OF_GENES),
                get_random(),
            ));
        }
    }

    pub fn sort_population(&mut self) {
        // https://rust-lang-nursery.github.io/rust-cookbook/algorithms/sorting.html
        self.population
            .sort_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap());
    }

    pub fn print_population(&mut self) {
        for p in 0..self.population.len() {
            println!("{}", self.population[p].fitness);
        }
    }

    pub fn item_crossover(&mut self, src_idx1: u32, src_idx2: u32) -> TIndivid {
        return crossover(
            self.population[src_idx1 as usize].to_owned(),
            self.population[src_idx2 as usize].to_owned(),
        );
    }

    pub fn check_limits(&mut self, dest_idx: u32) {
        self.population[dest_idx as usize].limit_values();
    }

    pub fn dump_to_file(&mut self, filepath: &str) {
        let temp_file = PathBuf::from(filepath); //party like it is 198x
                                                 // Open a file in write-only (ignoring errors).
                                                 // This creates the file if it does not exist (and empty the file if it exists).
        let mut file = File::create(temp_file).unwrap();
        let _res = writeln!(&mut file, "fitness,capacitance,resistance");
        for p in 0..self.population.len() {
            writeln!(
                &mut file,
                "{},{},{}",
                self.population[p].fitness,
                self.population[p].genes[0],
                self.population[p].genes[1]
            )
            .unwrap();
        }
        let mut _resultat = file.close(); //closing of files is handled different from most languages as a separate set of errors (see: https://docs.rs/close-file/latest/close_file/ )
    }
}

const ELITE_PART: f32 = 0.1;
const CROSSOVER_PART: f32 = 0.6;
const POPULATION_SIZE: u32 = 500;
const MAX_GENERATIONS: u32 = 20;
const MUTATION_INTENSITY: f32 = 0.5;
const MUTATION_RATE: f32 = 0.5;

pub(crate) fn main() {
    println!("Read CSV");
    println!("------------------------------------------------------------");
    let sim_data = read_csv("batterydata.csv");

    println!("Test of population_functions");
    println!("------------------------------------------------------------");
    let mut pop = TPopulation::new(MAX_GENERATIONS, MUTATION_RATE, MUTATION_INTENSITY);
    pop.make_population(POPULATION_SIZE); //population size

    let temp_file = PathBuf::from(r"fitness.csv"); //party like it is 198x
    let mut file = File::create(temp_file).unwrap();
    writeln!(&mut file, "generation, best_fitness, avg_fitness").unwrap();
    for current_generation in 0..pop.max_generation {
        println!(
            "- Generation {}----------------------------------------------------",
            current_generation
        );
        let mut avg_fitness: f32 = 0.0;
        let mut best_fitness: f32 = 1e9; //smaller is better (and always positive)
        for idx in 0..pop.population.len() {
            let fitness = pop.population[idx].calc_fitness(sim_data.to_owned(), "", false);
            //print!("< fit {},C{},R{}>",fitness,pop.population[idx].genes[0], pop.population[idx].genes[1]);
            if fitness < best_fitness {
                best_fitness = fitness;
            };
            avg_fitness = avg_fitness + fitness;
        }
        avg_fitness = avg_fitness / (pop.population.len() as f32);
        //println!(".");
        pop.sort_population();
        println!(
            "  current best: fitness={}, Capacitance = {} , Resistance={}",
            pop.population[0].fitness, pop.population[0].genes[0], pop.population[0].genes[1]
        );
        //pop.print_population();
        let elite_idx = (pop.population.len() as f32 * ELITE_PART) as u32;
        let crossover_idx = (pop.population.len() as f32 * CROSSOVER_PART) as u32;
        //let left_idx =pop.population.len() as u32;
        for idx in elite_idx..(pop.population.len() as u32) {
            //The part not elite is crossover
            let src_idx1 = (get_random() * elite_idx as f32) as u32;
            let src_idx2 = (get_random() * (pop.population.len() as f32)) as u32;
            pop.population[idx as usize] = pop.item_crossover(src_idx1, src_idx2);
            pop.population[idx as usize].limit_values();
        }
        for idx in crossover_idx..(pop.population.len() as u32) {
            //part of crossover is mutated
            pop.population[idx as usize].mutate(
                current_generation,
                pop.max_generation,
                pop.mutation_intensity,
                pop.mutation_rate,
            );
            pop.population[idx as usize].limit_values();
        }
        writeln!(
            &mut file,
            "{}, {}, {}",
            current_generation, best_fitness, avg_fitness
        )
        .unwrap();
    }

    //Store results
    pop.dump_to_file("last_generation.csv");
    println! {"Estimated capacitance {}, estimated resistance {}", pop.population[0].genes[0], pop.population[0].genes[1]};
    let mut _fitness = calc_voltage(
        "best_result.csv",
        true,
        sim_data.to_owned(),
        pop.population[0].genes[0],
        pop.population[0].genes[1],
    );
    println!("Done! (Caveat Emptor)");
}
