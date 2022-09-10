# Garry's mod ACF ballistic analysis with runge kutta
# trinket python: https://trinket.io/library/trinkets/e0aff07491

# ----- Space ----- #
# Set target distance
X = 10.0
Y = 1.0


# 3D world objects
box1  = box(   pos = vector(X/2, 0, 0), size   = vector(X, 0.2, 1), color = color.green)  # Ground
ball1 = sphere(pos = vector(0, 0, 0),   radius = 0.2,               color = color.yellow) # Shell
ball2 = sphere(pos = vector(X, Y, 0),   radius = 0.2,               color = color.red)    # Target


# ----- Setting ----- #
# Calc
CalcSpd   = 0.05  # Calculation speed
Tolerance = 0.25  # Tolerance of distance
LoopLimit = 50    # If iterations of solver loop is over than this value, will be stopping
Arc       = 0     # High, low fire (0: low, 1: high)

# Sys
G       = 600     # Gravity    (Source engine unit, 11.43m/s)
T       = 0.015   # Delta time

# Shell
D       = 40      # DragDiv    (ACF2: 40, ACF3: 80)
C       = 120     # Caliber    (mm unit)
M       = 0.01    # Mass
V       = 100     # Muzzle velocity

# PID gain
Kp      = 2       # Proportional
Ki      = 0.5     # Integral
Kd      = 0.1     # Differential
Kt      = 0.1     # Delta time


# ----- Constant ----- #
Cd = (pi*C**2)/(4000000*M*D)
T2 = T/2
T6 = T/6
Gv = vector(0, -G, 0)


# ----- Functions ----- #
# Return to radian
def toRad(th): return th*(pi/180)


# Return to degree
def toDeg(th): return th*(180/pi)


# Return vector2 length
def length(x, y): return sqrt(x**2 + y**2)


# Value min, max limiter
def clamp(n, min, max):
  if n < min: n = min
  if n > max: n = max
  return n


# Quadratic equation
def quadratic(sign):
  root = (V**4 - G*(G*X**2 + 2*Y*V**2))
  tan  = (V**2 + sign*sqrt(root))/(G*X)
  return toDeg(atan(tan))


# Runge Kutta 4
# Gravity accelation
def acc(k, v): return Gv - k*v


# Drag velocity
def vel(v, a, t): return v + a*t


def sv(t, v1, v2, v3, v4): return t*(v1 + 2*(v2 + v3) + v4)


def runge_kutta(th):
  pf = vector(0, 0, 0)
  vf = vector(V*cos(th), V*sin(th), 0)
  
  while True:
    rate(CalcSpd*1000)
    
    k  = length(vf.x, vf.y)*Cd
    v1 = vf
    a1 = acc(k, v1)
    v2 = vel(v1, a1, T2)
    a2 = acc(k, v2)
    v3 = vel(v1, a2, T2)
    a3 = acc(k, v3)
    v4 = vel(v1, a3, T)
    a4 = acc(k, v4)
    pf += sv(T6, v1, v2, v3, v4)
    vf += sv(T6, a1, a2, a3, a4)
    
    ball1.pos = pf
    limit     = (Y/X)*pf.x
    
    if pf.x > X | pf.y < limit:
      return pf


# Ballistic solver
def solve(arc, tolerance):
  if arc == 1: sign = 1
  else:        sign = -1
  
  # The first elements value
  ang  = quadratic(sign)
  th   = toRad(ang)
  loop = 0
  i    = pe = 0 # PID
  print("First angle:", ang)
  
  # Loop
  while abs(X - ball1.pos.x) + abs(Y - ball1.pos.y) >= tolerance:
    endPos = runge_kutta(th)
    
    # Error
    ex = X - endPos.x
    ey = Y - endPos.y
    er = clamp(ex + ey, -20, 20)
    
    # PID
    p   = er
    i   += p*Kt
    d   = (p - pe)/Kt
    pe  = er
    pid = p*Kp + i*Ki + d*Kd
    
    # Update elements
    ang  += pid*-sign
    th   = toRad(ang)
    loop += 1
    
    if ang > 90 | Y > 0 & ang < 0 | loop > LoopLimit:
      print("Solver is terminated")
      return th
  
  print("Optimal angle:", toDeg(th), "/ Loop:", loop)
  return th


# Solve
theta1    = solve(Arc, Tolerance)
ball1.pos = vector(0, 0, 0)

