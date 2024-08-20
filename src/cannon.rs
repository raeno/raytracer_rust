use crate::tuple::Tuple;
use std::fmt::Display;

type Point = Tuple<f64>;
type Vector = Tuple<f64>;

#[derive(Debug)]
pub struct Projectile {
    pub position: Point,
    pub velocity: Vector,
}

impl Display for Projectile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Position: ({}), Velocity: ({})",
            self.position, self.velocity
        )
    }
}

pub struct Environment {
    pub gravity: Vector,
    pub wind: Vector,
}

impl Projectile {
    pub fn on_surface(&self) -> bool {
        self.position.y >= 0.0
    }
}

pub fn projectile(position: Point, velocity: Vector) -> Projectile {
    Projectile { position, velocity }
}

pub fn environment(gravity: Vector, wind: Vector) -> Environment {
    Environment { gravity, wind }
}

pub fn tick(environment: &Environment, projectile: &Projectile) -> Projectile {
    let position = projectile.position + projectile.velocity;
    let velocity = projectile.velocity + environment.gravity + environment.wind;
    Projectile { position, velocity }
}
