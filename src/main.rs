pub mod tuple;
pub mod cannon;

use tuple::Normalize;

use crate::tuple::{point, vector};
use crate::cannon::{projectile, environment, tick};

fn main() {
    let mut bullet = projectile(point(0, 1, 0).into(), vector(1, 1, 0).normalize());
    let environment = environment(vector(0.0, -0.1, 0.0), vector(-0.01, 0.0, 0.0));
    while bullet.on_surface() {
        bullet = tick(&environment, &bullet);
        println!("Bullet position: {}", &bullet.position);
    }
    println!("Final bullet position: {}", bullet);
}

