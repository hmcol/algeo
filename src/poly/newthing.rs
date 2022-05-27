use std::env::Args;

struct Environment {
    vars: Vec<Variable>,
}

struct ZZ(i64);
struct QQ {
    num: ZZ,
    den: ZZ,
}

struct Func<Args, Ret> {
    call: dyn Fn(Args) -> Ret,
}

trait AddGroup {
    fn zero() -> Self;
    fn add(a: Self, b: Self) -> Self;
    fn neg(a: Self) -> Self;
}

trait MulGroup {
    fn one() -> Self;
    fn mul(a: Self, b: Self) -> Self;
    fn inv(a: Self) -> Self;
}

struct ArithmeticComputer;

trait Addition<A, B, C> {
    fn add(a: A, b: B) -> C;
}

impl Addition<Function, Function, Function> for ArithmeticComputer {
    fn add(a: Function, b: Function) -> Function {
        let mut args = a.args;
        for arg in b.args {
            if !args.contains(&arg) {
                args.push(arg);
            }
        }


        Function::Env(Box::new(
            EnvFn {
                args: todo!(),
                expr: todo!(),
            }
        ))
    }
}


struct Variable {
    name: String,
}

struct Constant {
    value: i64,
}

enum CoreFn {
    Add,
    Sub,
    Mul,
    Div,
    Neg,
    Inv,
    Pow,
    Abs,
}

enum Expr {
    Const(Constant),
    Var(Variable),
    Function(FnExpr),
}

struct EnvFn {
    args: Vec<Variable>,
    expr: Box<FnExpr>,
}

enum Function {
    Core(CoreFn),
    Env(EnvFn),
}

enum FnExpr {
    Var(Variable),
    Call(Function, Vec<FnExpr>),
}

impl FnExpr {
    fn eval(&self, vars: &Vec<Variable>) -> Expr {
        match self {
            FnExpr::Var(var) => Expr::Var(var.clone()),
            FnExpr::Call(func, args) => {
                let args = args.iter().map(|arg| arg.eval(vars)).collect();
                Expr::Function(FnExpr::Call(func.clone(), args))
            }
        }
    }
}

