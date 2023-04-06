use std::f32::consts::E;

use ndarray::prelude::*;
use ndarray::{arr1, arr2, s, Array, Array2, NewAxis};
use ndarray_linalg::*;
use rand_distr::num_traits::signum;

const z: f32 = 1.0; // Electron exchange number
const F: f32 = 96485.0; // Faraday constant
const R: f32 = 8.314; // Ideal gas constant
const E_ZERO_25C: f32 = 0.65;
const V_CELL: f32 = 0.0025 * 0.002; // Surface area times thickness
const V_CELL_DIFF:f32 = 3.4e-6;
pub struct ElectroChemModel {
    c_nominal: f32,          // Nominal total concentration
    pub c_cell: Array2<f32>, // Cell concentrations
    pub c_tank: Array2<f32>,     // Tank concentrations
    pub Voltage: f32,        // Stack voltage
    SOC: f32,                // State of Charge
    SOH: f32,                // State of Health
    V_cell: f32,             // Cell volume
    V_tank: f32,             // Tank volume
    N: f32,                  // Stack cell count
    S: f32,                  // Membrane surface area
    d: f32,                  // Membrane thickness
    D: Array2<f32>,          // Diffusion kinetics matrix
    E0: f32,                 // Formal potential
    ChargeOffset: f32,       //   
    DisChargeOffset: f32,    // 
    kp: f32,                 // Rate constant, anolyte
    kn: f32,                 // Rate constant, catholyte
    Rstack: f32,             // Ohmic stack resistance
    pub dt: f32,             // Sample time
    currentsigns: Array2<f32>,
}

impl ElectroChemModel {
    // Parameters to be estimated are diffusion_rate, rate_anolyte, rate_catholyte, stack_resistance
    pub fn new(
        nominal_anolyte: f32,
        nominal_catholyte: f32,
        diffusion_rate: f32,
        rate_anolyte: f32,
        rate_catholyte: f32,
        stack_resistance: f32,
        charge_offset:f32,
        discharge_offset:f32,
        sample_time: f32,
        isdiffusioncell: bool
    ) -> Self {
        if isdiffusioncell == false{
            Self {
                c_nominal: nominal_anolyte,
                c_cell: arr2(&[[nominal_anolyte], [nominal_catholyte], [0.0]]),
                c_tank: arr2(&[[nominal_anolyte], [nominal_catholyte], [0.0]]),
                SOC: 0.0,
                SOH: 1.0,
                V_cell: V_CELL, // 1.68e-4,
                V_tank: 50e-6,
                N: 1.0,
                S: 0.0025, //0.067,
                d: 33e-6,
                D: arr2(&[
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 2.0 * diffusion_rate],
                    [0.0, 0.0, -1.0 * diffusion_rate],
                ]),
                dt: sample_time,
                E0: E_ZERO_25C,
                ChargeOffset : charge_offset,
                DisChargeOffset : discharge_offset,
                kp: rate_anolyte,
                kn: rate_catholyte,
                Rstack: stack_resistance,
                currentsigns: arr2(&[[-1.0], [-1.0], [1.0]]),
                Voltage: 0.0,
            }
        }else{
            Self {
                c_nominal: nominal_anolyte,
                c_cell: arr2(&[[nominal_anolyte], [nominal_catholyte], [0.0]]),
                c_tank: arr2(&[[nominal_anolyte], [nominal_catholyte], [0.0]]),
                SOC: 0.0,
                SOH: 1.0,
                V_cell: V_CELL_DIFF, // 1.68e-4,
                V_tank: 50e-6,
                N: 1.0,
                S: 1e-4,
                d: 33e-6,
                D: arr2(&[
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 2.0 * diffusion_rate],
                    [0.0, 0.0, -1.0 * diffusion_rate],
                ]),
                dt: sample_time,
                E0: E_ZERO_25C,
                ChargeOffset : charge_offset,
                DisChargeOffset : discharge_offset,
                kp: rate_anolyte,
                kn: rate_catholyte,
                Rstack: stack_resistance,
                currentsigns: arr2(&[[-1.0], [-1.0], [1.0]]),
                Voltage: 0.0,
            }
            }
}
}

const MIN_CURRENT: f32 = 0.0001;
impl ElectroChemModel {
    pub fn TimeStep(
        &mut self,
        q: f32, //Flow [m^3/s]
        mut I: f32,
        dt: f32,
    ) -> f32 {
        self.dt = dt;
        let Q: Array2<f32> = arr2(&[[q, 0.0, 0.0], [0.0, q, 0.0], [0.0, 0.0, q]]);

        let mut currentpart: Array2<f32> = 1.0 / (z * F) * &self.currentsigns * I;

        let flowpart: Array2<f32> = Q.dot(&(&self.c_tank - &self.c_cell)) ;
        let diffpart: Array2<f32> = self.S / self.d * &self.D.dot(&self.c_cell);
        let dtcell = self.dt * (flowpart + diffpart + currentpart) / (self.V_cell);
        let dttank = self.dt * (self.N * Q.dot(&(&self.c_cell - &self.c_tank))) / self.V_tank;

        self.c_cell = &self.c_cell + dtcell;

        for ii in 0..self.c_cell.len(){
            if self.c_cell[[ii,0]] < 0.0{
                self.c_cell[[ii,0]] = 0.0;
            }
            if self.c_tank[[ii,0]] < 0.0{
                self.c_tank[[ii,0]] = 0.0;
            }
        }
        self.c_tank = &self.c_tank + dttank;

        // Butler-Volmer equations
        if (I < MIN_CURRENT) && (I > 0.0) {
            I = MIN_CURRENT;
        };
        if (I > -MIN_CURRENT) && (I < 0.0) {
            I = -MIN_CURRENT;
        };
        let jp: f32 = 1.0 / self.S
            * (F * self.kp * self.c_cell[[2, 0]].powf(0.5) * self.c_cell[[0, 0]].powf(0.5));
        let jn: f32 = 1.0 / self.S * (F * self.kn * self.c_cell[[1, 0]].powf(0.5)*1000.0.powf(0.5));

        let logtermp: f32 = 1.0 / (2.0 * jp * self.S) * I
            + ((1.0 / (2.0 * jp * self.S) * I).powf(2.0) + 1.0).sqrt();
        let logtermn: f32 = 1.0 / (2.0 * jn * self.S) * I
            + ((1.0 / (2.0 * jn * self.S) * I).powf(2.0) + 1.0).sqrt();

        let mut Vp: f32 = 2.0 * R * 333.15 / F * (logtermp + 1e-10).ln();
        if Vp.is_nan() || Vp.is_infinite() {
            Vp = 0.0;
        }

        let mut Vn: f32 = 2.0 * R * 333.15 / F * (logtermn + 1e-10).ln();

        if Vn.is_nan() || Vn.is_infinite() {
            Vn = 0.0;
        }

        let Vbv: f32 = Vp - Vn;

        let mut Vnernst:f32 = 0.0;

        
        Vnernst = R * 333.15 / F
            * (self.c_cell[[2, 0]] / (self.c_cell[[0, 0]])
                * (1000.0/(self.c_cell[[1, 0]] + 1.0e-11)))
                .ln();

        if Vnernst.is_nan() || Vnernst.is_infinite() {
            Vnernst = 0.0;
        }

//        self.Voltage = self.N * (self.E0 + Vbv + Vnernst + self.Rstack * self.S * I);
        let mut Voffset:f32 = 0.0;
        if I < 0.0{
            Voffset = self.DisChargeOffset;
        }else{
            Voffset = self.ChargeOffset;
        }
        self.Voltage = self.N * (self.E0 + Vbv + Vnernst + self.Rstack * I + Voffset);
        //println!("Stack voltage is {:?}",self.Voltage);
        //println!("Anolyte C1 concentration is {:?}",self.c_cell[[0,0]]);
        return self.Voltage;
    }

    pub fn TimeStep_DiffCell(
        &mut self,
        mut I: f32,
        dt: f32,
    ) -> f32 {
        self.dt = dt;

        let currentpart: Array2<f32> = 1.0 / (z * F) * &self.currentsigns * I;

        let diffpart: Array2<f32> = self.S / self.d * &self.D.dot(&self.c_cell);

        let dtcell = self.dt * (diffpart + currentpart) / (self.V_cell);
        self.c_cell = &self.c_cell + dtcell;
        for ii in 0..self.c_cell.len(){
            if self.c_cell[[ii,0]] < 0.0{
                self.c_cell[[ii,0]] = 0.0;
            }
            if self.c_tank[[ii,0]] < 0.0{
                self.c_tank[[ii,0]] = 0.0;
            }
        }
        self.c_tank.assign(&self.c_cell);

        // Butler-Volmer equations
        if (I < MIN_CURRENT) && (I > 0.0) {
            I = MIN_CURRENT;
        };
        if (I > -MIN_CURRENT) && (I < 0.0) {
            I = -MIN_CURRENT;
        };
        let jp: f32 = 1.0 / self.S * (F * self.kp * self.c_cell[[2, 0]].powf(0.5) * self.c_cell[[0, 0]].powf(0.5));
        let jn: f32 = 1.0 / self.S * (F * self.kn *1000.0.powf(0.5)*self.c_cell[[1, 0]].powf(0.5));

        let logtermp: f32 = 1.0 / (2.0 * jp * self.S) * I + ((1.0 / (2.0 * jp * self.S) * I).powf(2.0) + 1.0).sqrt();
        let logtermn: f32 = 1.0 / (2.0 * jn * self.S) * I + ((1.0 / (2.0 * jn * self.S) * I).powf(2.0) + 1.0).sqrt();

        let mut Vp: f32 = 2.0 * R * 333.15 / F * (logtermp).ln();
        if Vp.is_nan() || Vp.is_infinite() {
            Vp = 0.0;
        }

        let mut Vn: f32 = 2.0 * R * 333.15 / F * (logtermn).ln();

        if Vn.is_nan() || Vn.is_infinite() {
            Vn = 0.0;
        }

        let Vbv: f32 = Vp - Vn;

        let mut Vnernst: f32 = R * 333.15 / F * ((self.c_cell[[2, 0]]*1000.0)/(self.c_cell[[0, 0]]*self.c_cell[[1, 0]])).ln();

        if Vnernst.is_nan() || Vnernst.is_infinite() {
            Vnernst = 0.0;
        }

//        self.Voltage = self.N * (self.E0 + Vbv + Vnernst + self.Rstack * self.S * I);
        let mut Voffset:f32 = 0.0;
        if I < 0.0{
            Voffset = self.DisChargeOffset;
        }else{
            Voffset = self.ChargeOffset;
        }
        self.Voltage = self.N * (self.E0 + Vbv + Vnernst + self.Rstack * I + Voffset);
        //println!("Stack voltage is {:?}",self.Voltage);
        //println!("Anolyte C1 concentration is {:?}",self.c_cell[[0,0]]);
        return self.Voltage;
    }
}
