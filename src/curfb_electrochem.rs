use std::f32::consts::E;

use ndarray::{arr1, arr2, Array2, s, Array, NewAxis};
use rand_distr::num_traits::signum;
use ndarray::prelude::*;
use ndarray_linalg::*;

const z:f32 = 1.0; // Electron exchange number
const F:f32 = 96485.0; // Faraday constant
const R:f32 = 8.314; // Ideal gas constant

pub struct ElectroChemModel{
    c_nominal:f32, // Nominal total concentration
    c_cell:Array2<f32>, // Cell concentrations
    c_tank:Array2<f32>, // Tank concentrations
    pub Voltage:f32, // Stack voltage
    SOC:f32, // State of Charge
    SOH:f32, // State of Health
    V_cell:f32, // Cell volume
    V_tank:f32, // Tank volume
    N:f32, // Stack cell count
    S:f32, // Membrane surface area
    d:f32, // Membrane thickness
    D:Array2<f32>, // Diffusion kinetics matrix
    E0:f32, // Formal potential
    kp:f32, // Rate constant, anolyte
    kn:f32, // Rate constant, catholyte
    Rstack:f32, // Ohmic stack resistance
    dt:f32, // Sample time
    currentsigns:Array2<f32>
}

impl ElectroChemModel{ // Parameters to be estimated are diffusion_rate, rate_anolyte, rate_catholyte, stack_resistance 
    pub fn new(
        nominal_concentration:f32,diffusion_rate:f32,rate_anolyte:f32,
        rate_catholyte:f32,stack_resistance:f32, sample_time:f32
        ) -> Self{

        Self{
            c_nominal: nominal_concentration,
            c_cell: arr2(&[[1.0],[1.0],[0.0]])*nominal_concentration,
            c_tank: arr2(&[[1.0],[1.0],[0.0]])*nominal_concentration,
            SOC: 0.0,
            SOH: 1.0,
            V_cell: 9.5e-4,
            V_tank: 0.01,
            N: 1.0,
            S: 0.15,
            d: 1.27e-4,
            D: arr2(&[[0.0,0.0,0.0],[0.0,0.0,2.0*diffusion_rate],[0.0,0.0,-1.0*diffusion_rate],]),
            dt: sample_time,
            E0: 0.65,
            kp:rate_anolyte,
            kn:rate_catholyte,
            Rstack: stack_resistance,
            currentsigns: arr2(&[[-1.0],[-1.0],[1.0]]),
            Voltage: 0.0
        }
    }
}

impl ElectroChemModel{

    pub fn TimeStep(&mut self,q:f32,I:f32){
        let Q:Array2<f32> = arr2(&[
            [q,0.0,0.0],
            [0.0,q,0.0],
            [0.0,0.0,q],
            ]);

        let mut currentpart:Array2<f32> = 1.0/(z*F)*&self.currentsigns*I;

        let flowpart:Array2<f32> = Q.dot(&(&self.c_tank-&self.c_cell));

        let diffpart:Array2<f32> = self.S/self.d*&self.D.dot(&self.c_cell);

        let dtcell = self.dt*(flowpart+diffpart+currentpart)/(self.V_cell);
        let dttank = self.dt*(self.N*Q.dot(&(&self.c_cell-&self.c_tank)))/self.V_tank;
        self.c_cell = &self.c_cell + dtcell;
        self.c_tank = &self.c_tank + dttank;

        // Butler-Volmer equations
        let jp:f32 =  1.0/self.S*(F*self.kp*&self.c_cell[[2,0]].powf(0.5)*&self.c_cell[[0,0]].powf(0.5));
        let jn:f32 =  1.0/self.S*(F*self.kn*&self.c_cell[[1,0]].powf(0.5)); 
            
        let logtermp:f32 = 1.0/(2.0*jp*self.S)*I + ((1.0/(2.0*jp*self.S)).powf(2.0)+1.0).sqrt();
        let logtermn:f32 = 1.0/(2.0*jn*self.S)*I + ((1.0/(2.0*jn*self.S)).powf(2.0)+1.0).sqrt();

        let Vp:f32 = 2.0*R*323.15/F*logtermp.ln();
        let Vn:f32 = 2.0*R*323.15/F*logtermn.ln();
        let Vbv:f32 = Vp-Vn;

        let Vnernst:f32 = R*323.15/F*(self.c_cell[[2,0]]/(self.c_cell[[0,0]]+1.0e-11)*(1.0/(self.c_cell[[1,0]]+1.0e-11))).ln();

        self.Voltage = self.N*(self.E0+Vbv+Vnernst+self.Rstack*self.S*I);
        println!("Stack voltage is {:?}",self.Voltage);
        println!("Anolyte C1 concentration is {:?}",self.c_cell[[0,0]]); 

    }
}


// pub struct ElectroChemModel{
//     c_cell:Array2<f32>, // Cell concentrations
//     c_tank:Array2<f32>, // Tank concentrations
//     sample_time:f32, // Model sample time
//     currentsigns:Array2<f32>, // Current sign vector
//     A:Array2<f32> // Reaction kinetics matrix
// }

// impl ElectroChemModel{
//     pub fn new(c_cell_init:Array2<f32>,c_tank_init:Array2<f32>,sample_time:f32) -> Self{

//         Self{
//             c_cell: c_cell_init,
//             c_tank: c_tank_init,
//             sample_time: sample_time,
//             currentsigns: arr2(&[[-1.0],[-1.0],[0.0],[1.0]]),
//             A: arr2(&[[-k0,k0, 2.0*k_d, 0.0],[k0,-k0,0.0,0.0],[0.0,0.0,-(k_d-k1),k1],[0.0,0.0,k1,-k1]])
//         }

//     }
// }

// impl ElectroChemModel{

//     pub fn TimeStep(&mut self,qa:f32,qc:f32,I:f32){
//         let Q:Array2<f32> = arr2(&[
//             [qa,0.0,0.0,0.0],
//             [0.0,qc,0.0,0.0],
//             [0.0,0.0,qa,0.0],
//             [0.0,0.0,0.0,qc]
//             ]);
//         let mut currentpart:Array2<f32> = 1.0/(z*F)*&self.currentsigns*I;
//         let mut c0:f32 = 5000.0;
//         let sum:f32 = self.c_cell.iter().sum();
//         c0 = c0-sum;
//         if I > 0.01{
//             if (self.c_cell[[0,0]] < 0.001) && (self.c_cell[[1,0]] < 0.001){
//                 currentpart = currentpart*0.0;
//             }
//         }
//        else if I < -0.01  {
//             if (c0 < 0.001) && (self.c_cell[[3,0]] < 0.001){
//                 currentpart = currentpart*0.0;
//             } 
//         }   
//         else{
//             currentpart = currentpart*0.0;
//         }

//         let flowpart:Array2<f32> = Q.dot(&(&self.c_tank-&self.c_cell));
//         let diffpart:Array2<f32> = S/d*&self.A.dot(&self.c_cell);
//         let dtcell = self.sample_time*(flowpart+diffpart+currentpart)/(2.0*V_cell);
//         let dttank = self.sample_time*(N*Q.dot(&(&self.c_cell-&self.c_tank)))/V_tank;
//         self.c_cell = &self.c_cell + dtcell;
//         self.c_tank = &self.c_tank + dttank;
//     }
// }

// impl ElectroChemModel{
//     pub fn gettanks(&self) -> Array2<f32> {
//         let mut outvec:Array2<f32> = Array2::<f32>::zeros((self.c_tank.len(),1));
//         for ii in 0..self.c_tank.len()
//         {
//             outvec[[ii,0]] = self.c_tank[[ii,0]];
//         }
//         return outvec;
//     }

//     pub fn getcells(&self) -> Array2<f32> {
//         let mut outvec:Array2<f32> = Array2::<f32>::zeros((self.c_cell.len(),1));
//         for ii in 0..self.c_cell.len()
//         {
//             outvec[[ii,0]] = self.c_cell[[ii,0]];
//         }
//         return outvec;
//     }
// }
