import math
import copy
decimal = 11
def avg(a,b):
    return((a+b)/2)
def slope_intercept_linear_system(m1,b1,m2,b2):
    x = ((b2-b1)/(m1-m2))
    y = (m1*x+b1)
    return[x,y]
def standard_linear_system(a1,b1,c1,a2,b2,c2):
    x = (c1*b2-c2*b1)/(b1*a2-b2*a1)
    y = (c1*a2-c2*a1)/(a1*b2-a2*b1)
    return[x,y]
def two_point_line(x1,y1,x2,y2):
    m = (y1-y2)/(x1-x2)
    b = y1-m*x1
    return[m,b]
def point_slope_line(x,y,m):
    b = y-m*x
    return[m,b]
def point_intercept_line(x,y,b):
    m = (y-b)/x
    return[m,b]
def midpoint(x1,y1,x2,y2):
    x = (x1+x2)/2
    y = (y1+y2)/2
    return[x,y]
#point1 and point2 should be lists containing all the coordinates of the point. example input is ([1,2,3],[4,5,6]). both points should have the same number of coordinates
def distance(point1,point2):
    squared_sum = 0
    count = 0
    while(count<len(point1)):
        squared_sum += (point1[count]-point2[count])**2
        count += 1
    return(round(math.sqrt(squared_sum),decimal))
def perpendicular_bisector(x1,y1,x2,y2):
    x = (x1+x2)/2
    y = (y1+y2)/2
    if(x1 == x2):
        m = 0
    elif(y1 == y2):
        return(x)
    else:
        m = -1/((y1-y2)/(x1-x2))
    b = -1*m*x+y
    return [m,b]
def parallelogram_check(x1,y1,x2,y2,x3,y3,x4,y4):
    if(avg(x1,x2) == avg(x3,x4) and avg(y1,y2) == avg(y3,y4)):
        return (True,2)
    elif(avg(x1,x3) == avg(x2,x4) and avg(y1,y3) == avg(y2,y4)):
        return (True,3)
    elif(avg(x1,x4) == avg(x2,x3) and avg(y1,y4) == avg(y2,y3)):
        return (True,4)
    return False
def rectangle_check(x1,y1,x2,y2,x3,y3,x4,y4):
    if(parallelogram_check(x1,y1,x2,y2,x3,y3,x4,y4) == (True,2)):
        if(x1 == x4 or x1 == x3):
            if((y1 == y3 and x1 == x4) or (y1 == y4 and x1 == x3)):
                return True
            return False
        elif((((y1-y3)/(x1-x3))*((y1-y4))/(x1-x4)) == -1):
            return True
        return False
    elif(parallelogram_check(x1,y1,x2,y2,x3,y3,x4,y4) == (True,3)):
        if(x1 == x4 or x1 == x2):
            if((y1 == y2 and x1 == x4) or (y1 == y4 and x1 == x2)):
                return True
            return False
        elif((((y1-y2)/(x1-x2))*((y1-y4))/(x1-x4)) == -1):
            return True
        return False
    elif(parallelogram_check(x1,y1,x2,y2,x3,y3,x4,y4) == (True,4)):
        if(x1 == x2 or x1 == x3):
            if((y1 == y3 and x1 == x2) or (y1 == y2 and x1 == x3)):
                return True
            return False
        elif((((y1-y3)/(x1-x3))*((y1-y2))/(x1-x2)) == -1):
            return True
        return False
    return False
def rhombus_check(x1,y1,x2,y2,x3,y3,x4,y4):
    if(parallelogram_check(x1,y1,x2,y2,x3,y3,x4,y4) == (True,2)):
        if(distance([x1,y1],[x3,y3]) == distance([x2,y2],[x3,y3]) == distance([x2,y2],[x4,y4]) == distance([x1,y1],[x4,y4])):
            return True
        return False
    if(parallelogram_check(x1,y1,x2,y2,x3,y3,x4,y4) == (True,3)):
        if(distance([x1,y1],[x2,y2]) == distance([x2,y2],[x3,y3]) == distance([x3,y3],[x4,y4]) == distance([x1,y1],[x4,y4])):
            return True
        return False
    if(parallelogram_check(x1,y1,x2,y2,x3,y3,x4,y4) == (True,4)):
        if(distance([x1,y1],[x3,y3]) == distance([x4,y4],[x3,y3]) == distance([x2,y2],[x4,y4]) == distance([x1,y1],[x2,y2])):
            return True
        return False
    return False
def square_check(x1,y1,x2,y2,x3,y3,x4,y4):
    if(rectangle_check(x1,y1,x2,y2,x3,y3,x4,y4) == True and rhombus_check(x1,y1,x2,y2,x3,y3,x4,y4) == True):
        return True
    return False
def quadrilateral_categorizer(x1,y1,x2,y2,x3,y3,x4,y4):
    rectangle = False
    rhombus = False
    if(parallelogram_check(x1,y1,x2,y2,x3,y3,x4,y4) == False):
        return "nonparallelogram"
    elif(rectangle_check(x1,y1,x2,y2,x3,y3,x4,y4) == True and rhombus_check(x1,y1,x2,y2,x3,y3,x4,y4) == False):
        return "rectangle"
    elif(rectangle_check(x1,y1,x2,y2,x3,y3,x4,y4) == False and rhombus_check(x1,y1,x2,y2,x3,y3,x4,y4) == True):
        return "rhombus"
    elif(rectangle_check(x1,y1,x2,y2,x3,y3,x4,y4) == True == rhombus_check(x1,y1,x2,y2,x3,y3,x4,y4)):
        return "square"
    return "parallelogram"
def centroid(x1,y1,x2,y2,x3,y3):
    x = (x1+x2+x3)/3
    y = (y1+y2+y3)/3
    return[x,y]
def circumcenter(x1,y1,x2,y2,x3,y3):
    if(x1 == x2):
        m1 = 0
    elif(y1 == y2):
        m2 = -1/((y1-y3)/(x1-x3))
        b2 = -1*m2*(x1+x3)/2+(y1+y3)/2
        x = (x1+x2)/2
        y = m2*x+b2
        return[x,y]
    else:
        m1 = -1/((y1-y2)/(x1-x2))
    b1 = -1*m1*(x1+x2)/2+(y1+y2)/2
    if(x1 == x3):
        m2 = 0
    elif(y1 == y3):
        x = (x1+x3)/2
        y = m1*x+b1
        return[x,y]
    else:
        m2 = -1/((y1-y3)/(x1-x3))
    b2 = -1*m2*(x1+x3)/2+(y1+y3)/2
    return(slope_intercept_linear_system(m1,b1,m2,b2))
def orthocenter(x1,y1,x2,y2,x3,y3):
    if(x1 == x2):
        m1 = 0
    elif(y1 == y2):
        m2 = -1/((y1-y3)/(x1-x3))
        b2 = -1*m2*x2+y2
        x = x3
        y = m2*x+b2
        return[x,y]
    else:
        m1 = -1/((y1-y2)/(x1-x2))
    b1 = -1*m1*x3+y3
    if(x1 == x3):
        m2 = 0
    elif(y1 == y3):
        x = x2/2
        y = m1*x+b1
        return[x,y]
    else:
        m2 = -1/((y1-y3)/(x1-x3))
    b2 = -1*m2*x2+y2
    return(slope_intercept_linear_system(m1,b1,m2,b2))
def equilateral_check(x1,y1,x2,y2,x3,y3):
    if(distance([x1,y1],[x2,y2]) == distance([x1,y1],[x3,y3]) == distance([x3,y3],[x2,y2])):
        return True
    return False
def isosceles_check(x1,y1,x2,y2,x3,y3):
    if(distance([x1,y1],[x2,y2]) == distance([x1,y1],[x3,y3]) or distance([x1,y1],[x2,y2]) == distance([x3,y3],[x2,y2]) or distance([x1,y1],[x3,y3]) == distance([x3,y3],[x2,y2])):
        return [True]
    return [False]
def right_triangle_check(x1,y1,x2,y2,x3,y3):
    if(x1 == x2 or x1 == x3):
        if((y1 == y3 and x1 == x2) or (y1 == y2 and x1 == x3)):
           return True
    elif((((y1-y3)/(x1-x3))*((y1-y2))/(x1-x2)) == -1):
        return True
    if(x1 == x2 or x2 == x3):
        if((y2 == y3 and x1 == x2) or (y1 == y2 and x2 == x3)):
           return True
    elif((((y2-y3)/(x2-x3))*((y1-y2))/(x1-x2)) == -1):
        return True
    if(x3 == x2 or x1 == x3):
        if((y1 == y3 and x3 == x2) or (y3 == y2 and x1 == x3)):
           return True
    elif((((y1-y3)/(x1-x3))*((y3-y2))/(x3-x2)) == -1):
        return True
    return False
def triangle_categorizer(x1,y1,x2,y2,x3,y3):
    if(equilateral_check(x1,y1,x2,y2,x3,y3) == True):
        return "equilateral"
    elif(isosceles_check(x1,y1,x2,y2,x3,y3) == right_triangle_check(x1,y1,x2,y2,x3,y3) == True):
        return "right isosceles"
    elif(right_triangle_check(x1,y1,x2,y2,x3,y3) == True):
        return "right scalene"
    elif(isosceles_check(x1,y1,x2,y2,x3,y3) == True):
        return "non-right isosceles"
    return "non-right scalene"
def parabola_vertex(a,h,k):
    if((k<0 and a>0 ) or (k>0 and a<0)):
        root1 = ((2*a*h)-math.sqrt((2*a*h)**2-4*((a*h)**2+a*k)))/(2*a)
        root2 = ((2*a*h)+math.sqrt((2*a*h)**2-4*((a*h)**2+a*k)))/(2*a)
        return [h,k,[root1,root2]]
    elif(k == 0):
        return [h,k,[k]]
    else:
        return [h,k,[]]
def parabola_intercept(a,r,s):
    x = (r+s)/2
    y = a*(x-r)*(x-s)
    return[x,y,[r,s]]
def parabola_standard(a,b,c):
    k = c-(b**2/(4*a))
    h = -b/(2*a)
    if((k<0 and a>0 ) or (k>0 and a<0)):
        root1 = ((2*a*h)-math.sqrt((2*a*h)**2-4*((a*h)**2+a*k)))/(2*a)
        root2 = ((2*a*h)+math.sqrt((2*a*h)**2-4*((a*h)**2+a*k)))/(2*a)
        return [h,k,[root1,root2]]
    elif(k == 0):
        return[h,k,[h]]
    else:
        return [h,k,[]]
def three_point_parabola(x1,y1,x2,y2,x3,y3,format):
    a = (((y1-y2)/(x1-x2))-((y1-y3)/(x1-x3)))/(x2-x3)
    b = (y1-y2-a*(x1-x2)*(x1+x2))/(x1-x2)
    c = y1-a*x1**2-b*x1
    if(format == "s"):
        return [a,b,c]
    elif(format == "f"):
        points = parabola_standard(a,b,c)
        r = points[2][0]
        if(len(points[2]) == 1):
            return[a,r,r]
        else:
            s = points[2][1]
        return[a,r,s]
    elif(format == "v"):
        k = c-(b**2/(4*a))
        h = -b/(2*a)
        return [a,h,k]
def sidesideangle(a,b,A):
    conversion = input("radians(r) or degrees(d)?")
    print(conversion)
    if(conversion == "r"):
        conversion = 1
    elif(conversion == "d"):
        conversion = math.pi/180
    B = math.asin(math.sin(A*conversion)*b/a)/conversion
    C = math.pi/conversion-(B+A)
    c = a*math.sin(C*conversion)/math.sin(A*conversion)
    return[[a,A],[b,B],[c,C]]
def sideangleside(a,B,c):
    conversion = input("radians(r) or degrees(d)?")
    if(conversion == "r"):
        conversion = 1
    if(conversion == "d"):
        conversion = math.pi/180
    b = math.sqrt(a**2+c**2-2*a*c*math.cos(B*conversion))
    A = math.asin(math.sin(B*conversion)*a/b)/conversion
    C = math.pi/conversion-(A+B)
    return[[a,A],[b,B],[c,C]]
def angleangleside(A,B,a):
    conversion = input("radians(r) or degrees(d)?")
    if(conversion == "r"):
        conversion = 1
    if(conversion == "d"):
        conversion = math.pi/180
    C = math.pi/conversion-(A+B)
    b = a*math.sin(B*conversion)/math.sin(A*conversion)
    c = a*math.sin(C*conversion)/math.sin(A*conversion)
    return[[a,A],[b,B],[c,C]]
def anglesideangle(A,b,C):
    conversion = input("radians(r) or degrees(d)?")
    if(conversion == "r"):
        conversion = 1
    if(conversion == "d"):
        conversion = math.pi/180
    B = math.pi/conversion-(A+C)
    a = b*math.sin(A*conversion)/math.sin(B*conversion)
    c = b*math.sin(C*conversion)/math.sin(B*conversion)
    return[[a,A],[b,B],[c,C]]
def sidesideside(a,b,c):
    conversion = input("radians(r) or degrees(d)?")
    if(conversion == "r"):
        conversion = 1
    if(conversion == "d"):
        conversion = math.pi/180
    A = math.acos((b**2+c**2-a**2)/(2*b*c))/conversion
    B = math.acos((a**2+c**2-b**2)/(2*a*c))/conversion
    C = math.pi/conversion-(A+B)
    return[[a,A],[b,B],[c,C]]
def triangle_area(x1,y1,x2,y2,x3,y3):
    if(x1 == x2):
        height = math.sqrt((x3-x1)**2)
        base = distance([x1,y1],[x2,y2])
        area = base*height/2
        return [area]
    elif(y1 == y2):
        height = math.sqrt((y3-y1)**2)
        base = distance([x1,y1],[x2,y2])
        area = base*height/2
        return [area]
    else:
        m1 = -1/((y1-y2)/(x1-x2))
        b1 = -1*m1*x3+y3
    m2 = (y2-y1)/(x2-x1)
    b2 = (m2*x2-y2)*-1
    x = ((b2-b1)/(m1-m2))
    y = (m1*x+b1)
    base = distance([x1,y1],[x2,y2])
    height = distance([x3,y3],[x,y])
    area = round(base*height/2,decimal)
    return area
def four_point_cubic(x1,y1,x2,y2,x3,y3,x4,y4):
    a = (((y1-y2)/(x1-x2)-(y1-y3)/(x1-x3))/(x2-x3)-((y1-y2)/(x1-x2)-(y1-y4)/(x1-x4))/(x2-x4))/(x3-x4)
    b = (((y1-y2+a*(x2**3-x1**3))/(x1-x2))-((y1-y3+a*(x3**3-x1**3))/(x1-x3)))/(x2-x3)
    c = (y1-y2+a*(x2**3-x1**3)+b*(x2**2-x1**2))/(x1-x2)
    d = y1-a*x1**3-b*x1**2-c*x1
    return[a,b,c,d]
#will take two points (x1,y1) and (x2,y2) and return a and b for the function f(x)=a(b^x). NOTE: depending on how I feel during and after taking calculus, I might make it output a and k for the function f(x) = a(e^(kx)). If that's what you need, no worries, k = ln(b) and b = e^k so just do that to switch between them
def two_point_exponential(x1,y1,x2,y2):
    b = math.pow(math.e,math.log(y1/y2)/(x1-x2))
    a = y1/math.pow(b,x1)
    return[a,b]
#returns the a, b, and c values of a parabola that passes through points (x1,y1) and (x2,y2) with slope m at (x2,y2)
def two_point_slope_parabola(x1,y1,x2,y2,m):
    x = x2-x1
    x3 = x2+x
    y3 = m*x3+point_slope_line(x1,y1,m)[1]
    return(three_point_parabola(x1,y1,x2,y2,x3,y3,"s"))
def quadratic_system(a1,b1,c1,a2,b2,c2):
    a = a1-a2
    b = b1-b2
    c = c1-c2
    x = []
    y = []
    i = 0
    while(i<len(parabola_standard(a,b,c)[2])):
        x.append(parabola_standard(a,b,c)[2][i])
        y.append(a1*pow(x[-1],2)+b1*x[-1]+c1)
        i += 1
    return[x,y]
def three_point_circle(x1,y1,x2,y2,x3,y3):
    center = circumcenter(x1,y1,x2,y2,x3,y3)
    dist = distance([center[0],center[1]],[x1,y1])[0]
    return[center[0],center[1],dist]
def linear_regression(points):
    meanx = 0
    meanx2 = 0
    meanxy = 0
    meany = 0
    meany2 = 0
    squared_line_y_error = 0
    squared_mean_y_error = 0
    for i in points:
        meanx += i[0]
        meanx2 += i[0]**2
        meanxy += i[0]*i[1]
        meany += i[1]
        meany2 += i[1]**2
    meanx = meanx/len(points)
    meanx2 = meanx2/len(points)
    meanxy = meanxy/len(points)
    meany = meany/len(points)
    meany2 = meany2/len(points)
    x1 = meanx
    y1 = meany
    x2 = meanx2/meanx
    y2 = meanxy/meanx
    m = two_point_line(x1,y1,x2,y2)[0]
    b = two_point_line(x1,y1,x2,y2)[1]
    for j in points:
        squared_line_y_error += (j[1]-(m*j[0]+b))**2
        squared_mean_y_error += (j[1]-meany)**2
    r_squared = 1-squared_line_y_error/squared_mean_y_error
    return[[m,b],[r_squared]]
def quadratic_regression(points):
    meanx = 0
    meanx2 = 0
    meanx3 = 0
    meanx4 = 0
    meanxy = 0
    meanx2y = 0
    meany = 0
    meany2 = 0
    for i in points:
        meanx += i[0]
        meanx2 += i[0]**2
        meanx3 += i[0]**3
        meanx4 += i[0]**4
        meanxy += i[0]*i[1]
        meanx2y += i[0]**2*i[1]
        meany += i[1]
        meany2 += i[1]**2
    meanx = meanx/len(points)
    meanx2 = meanx2/len(points)
    meanx3 = meanx3/len(points)
    meanx4 = meanx4/len(points)
    meanxy = meanxy/len(points)
    meanx2y = meanx2y/len(points)
    meany = meany/len(points)
    meany2 = meany2/len(points)
    a1 = meanx2-meanx3/meanx
    a2 = meanx2-meanx4/meanx2
    b1 = meanx-meanx2/meanx
    b2 = meanx-meanx3/meanx2
    c1 = meanxy/meanx-meany
    c2 = meanx2y/meanx2-meany
    a = standard_linear_system(a1,b1,c1,a2,b2,c2)[0]
    b = standard_linear_system(a1,b1,c1,a2,b2,c2)[1]
    c = meany-a*meanx2-b*meanx
    R = meany2-2*meanx2y*a-2*meanxy*b-2*meany*c+a**2*meanx4+2*a*b*meanx3+2*a*c*meanx2+2*b*c*meanx+b**2*meanx2+c**2
    return[[a,b,c],[R]]
def exponential_regression(points):
    points1 = points.copy()
    if(len(points1[0]) == 1):
        year = 0
        while(year<len(points)):
            points1[year].insert(0,year)
            year = year+1
    i = 0
    #goes through, takes the natural log of each y value of each point in the array(the y value should be the 2nd value for each point, [x,y] NOT [y,x]
    while(i<len(points1)):
        points1[i][1] = math.log(abs(points1[i][1]),math.e)*abs(points1[i][1])/points1[i][1]
        i = i+1
    #then takes a linear regression of the new array
    lnline = linear_regression(points1)
    b = pow(math.e,lnline[0][0])
    a = pow(math.e,lnline[0][1])
    R = lnline[1][0]
    return[[a,b],[R]]
#like the exponential regression, but first takes the points and weights them algebraically by copying points based on position(0 will appear once, 1 twice, 6 seven times, and such) to weight the regression to more recent events. It is absolutely critical that you order the array chronologically, with more recent events at the end, so that it getss weighted properly
def weighted_exponential_regression(points):
    points1 = points.copy()
    year = 0
    while(year<len(points)):
        points1[year].insert(0,year)
        year = year+1
    i = 0
    while(i<len(points)):
        j = 0
        while(j<i):
            points1.append(points[i].copy())
            j = j+1
        i = i+1
    return(exponential_regression(points1))
def discounted_cashflow(earnings,expected_increase,discount_rate,years,historical_PE):
    future_value = earnings*pow(expected_increase,years)*historical_PE
    current_value = future_value/pow(discount_rate,years)
    return[current_value]
def value(yearly_equity,historical_PE,earnings,discount_rate,years):
    #the exponential regression function changes the array, so i make 2 arrays to use later and not interfere with each other
    weighted_yearly_equity = copy.deepcopy(yearly_equity)
    unweighted_yearly_equity = copy.deepcopy(yearly_equity)
    #so it takes both a weighted and unweighted exponential regression, the purposes of which will be explained shortly
    unweighted_equation = exponential_regression(yearly_equity)
    unweighted_expected_increase = unweighted_equation[0][1]
    weighted_equation = weighted_exponential_regression(weighted_yearly_equity)
    weighted_expected_increase = weighted_equation[0][1]
    #so i had this idea that, to factor in growth slowdown as companies get larger, you would take the unweighted average(i.e. 1.35 for 35%) and the weighted average(i.e. 1.3 for 30%), if it was smaller(which would indicate a growth slowdown, because the frontweightedness would make the data closer to currentish numbers, and so a lower weighted average would mean that the past few years have been slower than average, ergo its slowing down), you would take the difference between them and subtract that from the weighted one(i.e. 1.35-1.3 = 0.05, 1.3-0.05 = 1.25 for 25%) and that would account for slowdown slightly
    if(weighted_expected_increase<unweighted_expected_increase):
        expected_increase = weighted_expected_increase-(unweighted_expected_increase-weighted_expected_increase)
    else:
        expected_increase = weighted_expected_increase
    value = discounted_cashflow(earnings,expected_increase,discount_rate,years,historical_PE)
    return[value,expected_increase]
def mortgage_calc(principal,yearly_interest_rate,year_amount,periods_per_year):
    interest_rate = pow((1+yearly_interest_rate),1/periods_per_year)-1
    period_amount = year_amount*periods_per_year
    periodic_payment = (principal*pow((1+interest_rate),period_amount)*(-interest_rate))/(1-pow((1+interest_rate),(period_amount-1)))
    return[periodic_payment]
def stock_vs_house(current_cash,annual_stock_returns,capital_gains_tax_rate,annual_house_growth,annual_rental_income,annual_mortgage_interest_rate,mortgage_length_years,housing_expenses_without_mortgage,house_price,house_trade_expense,periods_per_year):
    stock_wealth_pretax = current_cash*pow(annual_stock_returns+1,mortgage_length_years)
    stock_tax = (stock_wealth_pretax-current_cash)*capital_gains_tax_rate
    stock_wealth = stock_wealth_pretax-stock_tax
    mortgage_principal = house_price+house_trade_expense-current_cash
    yearly_mortgage_payment = mortgage_calc(mortgage_principal,annual_mortgage_interest_rate,mortgage_length_years,periods_per_year)[0]*periods_per_year
    house_wealth = mortgage_length_years*(annual_rental_income-housing_expenses_without_mortgage-yearly_mortgage_payment)+house_price*pow(annual_house_growth,mortgage_length_years)-house_trade_expense
    return[stock_wealth,house_wealth]
def vector_addition(magnitude1,degrees_from_east1,magnitude2,degrees_from_east2):
    #this thing takes two vectors and adds them together. note that it will only return angles towards the east, so do actual math to see whether the resulting angle will be east or west(dw, you can probably do this with rough estimation) if it is on the wrong side, just add 180 degrees to the angle and itll be right
    conversion = math.pi/180
    x = magnitude1*math.cos(degrees_from_east1*conversion)+magnitude2*math.cos(degrees_from_east2*conversion)
    y = magnitude1*math.sin(degrees_from_east1*conversion)+magnitude2*math.sin(degrees_from_east2*conversion)
    magnitude = math.sqrt(pow(x,2)+pow(y,2))
    direction = math.atan(y/x)/conversion
    return[magnitude, direction]
def kinematic_equations(vi,vf,a,d,t):
    #stands for initial velocity, final velocity, acceleration(usually -9.8), displacement, change in time. put 3 of the values in in their rightful places, then fill in the other 2 values with None the syntax is important here so spell it None. also for projectile motion if you just have 1 velocity, acceleration, and displacement, it will only return values for the point BEFORE the parabola's apex, so keep that in mind and maybe dont use this for that
    if(vi == None and vf == None):
        vi1 = (d-0.5*a*pow(t,2))/t
        vf1 = (d+0.5*a*pow(t,2))/t
        a1 = a
        d1 = d
        t1 = t
    if(vi == None and a == None):
        vi1 = 2*d/t-vf
        vf1 = vf
        a1 = (vf-vi1)/t
        d1 = d
        t1 = t
    if(vi == None and d == None):
        vi1 = vf-a*t
        vf1 = vf
        a1 = a
        d1 = (vi1+vf)*t/2
        t1 = t
    if(vi == None and t == None):
        #can only be positive, so might be wrong sometimes
        vi1 = math.sqrt(pow(vf,2)-2*a*d)
        vf1 = vf
        a1 = a
        d1 = d
        t1 = (vf-vi1)/a
    if(vf == None and a == None):
        vi1 = vi
        vf1 = 2*d/t-vi
        a1 = (vf1-vi)/t
        d1 = d
        t1 = t
    if(vf == None and d == None):
        vi1 = vi
        vf1 = vi+a*t
        a1 = a
        d1 = (vi+vf1)*t/2
        t1 = t
    if(vf == None and t == None):
        vi1 = vi
        #can only be positive, so might be wrong sometimes
        vf1 = math.sqrt(pow(vi,2)+2*a*d)
        a1 = a
        d1 = d
        t1 = (vf1-vi)/a
    if(a == None and d == None):
        vi1 = vi
        vf1 = vf
        a1 = (vf-vi)/t
        d1 = (vi+vf)*t/2
        t1 = t
    if(a == None and t == None):
        vi1 = vi
        vf1 = vf
        a1 = (pow(vf,2)-pow(vi,2))/(2*d)
        d1 = d
        t1 = (vf-vi)/a1
    if(d == None and t == None):
        vi1 = vi
        vf1 = vf
        a1 = a
        d1 = (pow(vf,2)-pow(vi,2))/(2*a)
        t1 = (vf-vi)/a
    return[vi1,vf1,a1,d1,t1]
def retirement_number(years_left, expected_return,expected_inflation_plus_fees,yearly_contribution,retirement_expenses):
    #both expected_return and expected_inflation_plus_fees should be multipliers, like 1.06 and 0.96
    net_return = expected_return*expected_inflation_plus_fees
    year = 0
    retirement_cash = 0
    years_of_retirement = 0
    yearly_pension = (607.46+679.16)*12
    while(year<years_left):
        retirement_cash += yearly_contribution
        retirement_cash = retirement_cash*net_return
        year += 1
    cash_left = retirement_cash
    input_cash = yearly_contribution*years_left
    while(cash_left>0 and years_of_retirement<100):
        tax = (retirement_expenses-yearly_pension)*0.5
        if(cash_left>input_cash and years_of_retirement<100):
            cash_left = cash_left*net_return-retirement_expenses+yearly_pension-tax
        else:
            cash_left = cash_left*net_return-retirement_expenses+yearly_pension
        years_of_retirement += 1
    return[retirement_cash,years_of_retirement]
#The following three functions are based off of Invested by Phil & Danielle Town, for their pricing companies section in the middle of the book or something
#First one, tenhead pricing, it's when you take a company's owner earnings(net operating cash on cash flow statement, then unfactor tax and factor in maintenance expenditures(should be an add and subtract)) and multiply it by 10
#all variables should be able to be found on a 10-K, except maintenance_capital_expenditures. That has to be estimated, based on what YOU think the business has to spend on capital to maintain its current operations, as in how much it needs to spend to have 0% growth. this might be able to be found in footnotes of a 10-K, or it might be labelled under 'other' property & equipment expenses
#it takes those and returns a price that's 10*owner earnings, so that even with no company growth, you get a 10% return per year. EZ!!!
#tax should be positive, expenditure negative
def ten_cap_price(net_operating_cash,income_tax,maintenance_capital_expenditures):
    owner_earnings = net_operating_cash+income_tax+maintenance_capital_expenditures
    ten_cap_price = owner_earnings*10
    return[ten_cap_price]
#This one is kinda like tencap, but it takes into account compounding growth and spits out the price you would need to be paid back in 8 years of free cash flow. so like you would pay X dollars for company, in which you would expect to get all that money back in 8 years through the free cash flows added together
#inputs are all on financial statements, except maintenance & growth capital expenditures and expected growth rate
#maintenance & growth capital expenditures are similar to the one for tencap, but this time it's the expenditures a business needs to maintain and grow. Should be larger than just maintenance expenses. the idea is this+property_equipment purchases is going to be all the extra expenses the company needs to grow at your estimated rate, so it won't be going to free cash flow. both numbers should be in investing section of cash flow statement
#growth rate is the difficult one. this one is an estimate of how much you think the company's gonna grow year over year(i.e. 1.03 for 3% yearly growth) and it's mad inaccurate. have this estimate on the conservative side
#so what this does is it calculates free cash flow as operating cash(positive)+property_equipment_purchases(negative)+maintenance&growth capital expenditures(negative) and then does a lil geometric series sum for the free cash flows of the next 8 years compounded, and boom. das it
def eight_year_payback_price(expected_growth_rate,net_operating_cash,maintenance_and_growth_capital_expenditures):
    free_cash_flow = net_operating_cash+maintenance_and_growth_capital_expenditures
    eight_year_payback_price = free_cash_flow*((1-pow(expected_growth_rate,9))/(1-expected_growth_rate))-free_cash_flow
    return[eight_year_payback_price]
#Step 1: take earnings per share(EPS) and the expected_growth_rate(i.e. 1.04 for 4% a year) and calculate the earnings per share in ten years based off of that
#Step 2: multiply earnings by price to earnings ratio(P/E). you get P/E by getting twice the expected growth rate*2(i.e. 4%*2 = a P/E of 8) and the highest P/E over the last ten years, and picking the lower of the two
#Step 3: divide future price by 4, so we get 4 times our money in 10 years, an annualized return of about 15%
#Step 4: divide by 2 for that extra margin of safety
def discounted_cashflow_MOS(EPS,expected_growth_rate,highest_PE_over_last_10_years):
    future_earnings = EPS*pow(expected_growth_rate,10)
    if(highest_PE_over_last_10_years < 200*(expected_growth_rate-1)):
        PE = highest_PE_over_last_10_years
    else:
        PE = 200*(expected_growth_rate-1)
    future_price = future_earnings*PE
    current_price = future_price/4
    price_with_MOS = current_price/2
    return[price_with_MOS]
#pretty simple, just runs the last 3 programs together in one ez package, except you put net_income instead of EPS, so you get a full company valuation instead of a per share price to make comparison easier
def stock_price_calculator_three_ways(expected_growth_rate,net_operating_cash,income_tax,maintenance_and_growth_capital_expenditures,net_income,highest_PE_over_last_10_years):
    tencap_price = ten_cap_price(net_operating_cash,income_tax,maintenance_and_growth_capital_expenditures)
    eightyear_payback_price = eight_year_payback_price(expected_growth_rate,net_operating_cash,maintenance_and_growth_capital_expenditures)
    discounted_cash_flow = discounted_cashflow_MOS(net_income,expected_growth_rate,highest_PE_over_last_10_years)
    return[tencap_price,eightyear_payback_price,discounted_cash_flow]
#like the value function, but instead of valuing the company for a 15% return, it takes the current market cap and estimates what return you would get if you bought now and the future looks like the past
#note. don't use value and estimated annualized return on the same array, because the function changes the array. i might fix that
def estimated_annualized_return(yearly_equity,historical_PE,earnings,current_price,years):
    #the exponential regression function changes the array, so i make 2 arrays to use later and not interfere with each other
    weighted_yearly_equity = copy.deepcopy(yearly_equity)
    unweighted_yearly_equity = copy.deepcopy(yearly_equity)
    #so it takes both a weighted and unweighted exponential regression, the purposes of which will be explained shortly
    unweighted_equation = exponential_regression(yearly_equity)
    unweighted_expected_increase = unweighted_equation[0][1]
    weighted_equation = weighted_exponential_regression(weighted_yearly_equity)
    weighted_expected_increase = weighted_equation[0][1]
    #so i had this idea that, to factor in growth slowdown as companies get larger, you would take the unweighted average(i.e. 1.35 for 35%) and the weighted average(i.e. 1.3 for 30%), if it was smaller(which would indicate a growth slowdown, because the frontweightedness would make the data closer to currentish numbers, and so a lower weighted average would mean that the past few years have been slower than average, ergo its slowing down), you would take the difference between them and subtract that from the weighted one(i.e. 1.35-1.3 = 0.05, 1.3-0.05 = 1.25 for 25%) and that would account for slowdown slightly
    if(weighted_expected_increase<unweighted_expected_increase):
        expected_increase = weighted_expected_increase-(unweighted_expected_increase-weighted_expected_increase)
    else:
        expected_increase = weighted_expected_increase
    future_price = earnings*pow(expected_increase,years)*historical_PE
    increase = future_price/current_price
    estimated_annualized_return_rate = pow(increase,(1/years))
    return[estimated_annualized_return_rate,expected_increase]
#im such a big brain. so what you do here is input a list of the daily prices of a stock you wanna short term trade, i.e. [128.4,127.6,...,129.3] for like several years. then you input how much money you have to buy this stock, and how much it costs you to make each trade. finally, it asks for the maximum time frame moving average(MA) you want to test for(i.e. 21 for 21 days as the largest moving average to be tested) then, it simulates you buying and selling that stock with two moving averages as an indicator for the duration of the data. for example, imagine you want to buy when the 3-day MA is higher than the 7-day MA. this will simulate you buying and selling based on that indicator, and for all the other combinations of MA lengths up to the maximum time frame, and tell you which ORDERED pair of moving averages WOULD HAVE given you the best return, IF you had used it for the period of the stock which you inputted into the program(this is why i recommend a very long timeframe, like the past 15 years of daily close price for an s&p 500 index), considering the fact that you had to pay the fee every single time you bought or sold
#extra things: it assumes that each time you get a 'buy' or 'sell' signal, you would buy or sell all of the stock you can. 100%.
#also extra: it will generate a massive list of the results of all possible strategies used, but it will only print the best one and the "default" one, which is buying on day one and selling on the last day.
def best_simple_moving_average(daily_prices,starting_capital,fee,maximum_MA_time_frame):
    #test phase 1 is where MAlength2>MAlength1, phase 2 is the opposite, so i test a 3/7 day indicator, and it's opposite. i need to do both phases rather than just checking the lowest indicators, because the fee changes things. phase 3 means it's over
    test_phase = 1
    #when cash is true, you can buy but not sell, and if cash is false, you can sell but not buy. cash is to simulate whether your money is in cash or in stocks at the moment
    cash = True
    #setting the initial moving averages as 1 day and 2 day, along with their associated sentinel values
    MAlength1 = 1
    MAcounter1 = 0
    MAlength2 = 2
    MAcounter2 = 0
    #initiates ending capital as starting capital(the capital will change with each trade) and declaring both moving average prices
    ending_capital = starting_capital
    MA1 = 0
    MA2 = 0
    #initiates the list where the results from each tests will go. the list will eventually fill up with entries in the form of [[MAlength1,MAlength2],ending_capital] and be sorted by ending_capital so you can see the best indicator
    results = []
    while(test_phase != 3):
        #starts test at the earliest possible point in the list so that both moving averages will be fully utilized(i.e, if it's a 16/4 test, we start at [15] so that the 16 day MA has 16 data points to calculate
        if(MAlength2>MAlength1):
            current_day = MAlength2
        else:
            current_day = MAlength1
        MA1 = 0
        MA2 = 0
        MAcounter1 = 1
        MAcounter2 = 1
        while(MAcounter1<=MAlength1):
            MA1 += daily_prices[current_day-MAcounter1]
            MAcounter1 += 1
        while(MAcounter2<=MAlength2):
            MA2 += daily_prices[current_day-MAcounter2]
            MAcounter2 += 1
        #this loop runs through all the days in the daily prices list
        while(current_day<len(daily_prices)):
            #these two loops set the moving average values. it adds the price to the moving average for the correct amount of days, then divides by the amount of days
            #check whether to buy or sell
            MA1 += (daily_prices[current_day]-daily_prices[current_day-MAlength1])
            MA2 += (daily_prices[current_day]-daily_prices[current_day-MAlength2])
            if((MA1/MAlength1)>(MA2/MAlength2) and cash == True):
                cash = False
                ending_capital = (ending_capital-fee)/daily_prices[current_day]
            elif((MA1/MAlength1)<(MA2/MAlength2) and cash == False):
                cash = True
                ending_capital = ending_capital*daily_prices[current_day]-fee
            current_day += 1
        if(cash == False):
            cash = True
            ending_capital = ending_capital*daily_prices[current_day-1]-fee
        #this part adds the result to the end of the list, then changes the length of the moving averages being tested to check aaaaaaallllll the combos
        results.append([ending_capital,[MAlength1,MAlength2]])
        ending_capital = starting_capital
        if(test_phase == 1):
            MAlength2 += 1
            if(MAlength2>maximum_MA_time_frame):
                print(len(results))
                MAlength1 += 1
                MAlength2 = MAlength1 + 1
            if(MAlength1 == maximum_MA_time_frame):
                test_phase = 2
                print("halftime")
                MAlength1 = 2
                MAlength2 = 1
        elif(test_phase == 2):
            MAlength1 += 1
            if(MAlength1>maximum_MA_time_frame):
                print(len(results))
                MAlength2 += 1
                MAlength1 = MAlength2 + 1
            if(MAlength2 == maximum_MA_time_frame):
                test_phase = 3
    results.append([starting_capital/daily_prices[0]*daily_prices[len(daily_prices)-1],"Default"])
    sorted_results = sorted(results)
    return[sorted_results]
#ok this is kinda the same idea as best_simple_moving_average but using exponential and weighted moving averages as well. it tests in 2 phases, each with 6 tests, 3 intra-strategy tests and 3 inter-strategy tests done both forward and in reverse: simple moving average 1 against simple moving average 2, exponential moving average 1 against exponential moving average 2, weighted moving average 1 against weighted moving average 2 (those were the intras) and simple moving average 1 against exponential moving average 1, exponential moving average 1 against weighted moving average 1, and simple moving average 1 against weighted moving average 1 (those were the inters)
def best_moving_average(daily_prices,starting_capital,fee,maximum_MA_time_frame):
    #test phase 1 is where MAlength2>MAlength1, phase 2 is the opposite, so i test a 3/7 day indicator, and it's opposite. i need to do both phases rather than just checking the lowest indicators, because the fee changes things. phase 7 means it's over
    test_phase = 1
    #when cash is true, you can buy but not sell, and if cash is false, you can sell but not buy. cash is to simulate whether your money is in cash or in stocks at the moment
    #IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! the order of the cash and ending_capital lists goes simple vs simple, exponential vs exponential, weighted vs weighted, simple vs exponential, exponential vs weighted, and simple vs weighted
    cash = [True,True,True,True,True,True]
    #setting the initial moving averages as 1 day and 2 day, along with their associated sentinel values
    #NOTE: the exponential moving average will use these length values in a strange way. the exponential moving average formula ema(i)=ema(i-1)*(1-m)+A[i]*m, where m is the multiplicitave smoothing constant. Normal people use m=2/(time frame+1), so for a 10 day moving average m=1/11 = about 9%. This function, for simplicity, will not operate quite like this and will use m=1/(MAlength)
    MAlength1 = 1
    MAcounter1 = 0
    MAlength2 = 2
    MAcounter2 = 0
    #initiates ending capitals as starting capital(the capitals will change with each trade) and declaring all 6 moving average prices
    ending_capital = [starting_capital,starting_capital,starting_capital,starting_capital,starting_capital,starting_capital]
    SMA1 = 0
    SMA2 = 0
    EMA1 = 0
    EMA2 = 0
    WMA1 = 0
    WMA2 = 0
    #initiates the list where the results from each tests will go. the list will eventually fill up with entries in the form of [ending_capital,[MAlength1,MAlength2],test_type] and be sorted by ending_capital so you can see the best indicator
    results = []
    while(test_phase != 3):
        #starts test at the earliest possible point in the list so that all moving averages will be fully utilized(i.e, if it's a 16/4 test, we start at [15] so that the 16 day MA has 16 data points to calculate. This does not apply to the exponential moving average, which will always start at day 1(although this change of starting time should be insignificant for the vast majority of data points because the earlier data points become a smaller and smaller part of the average over time
        if(MAlength2>MAlength1):
            current_day = MAlength2
        else:
            current_day = MAlength1
        SMA1 = 0
        SMA2 = 0
        EMA1 = daily_prices[0]
        EMAM1 = 1/MAlength1
        EMA2 = daily_prices[0]
        EMAM2 = 1/MAlength2
        WMA1 = 0
        WMA2 = 0
        MAcounter1 = 1
        MAcounter2 = 1
        #these loops set up the initial values of all 6 averages on the current day
        while(MAcounter1<=MAlength1):
            SMA1 += daily_prices[current_day-MAcounter1]
            EMA1 = EMA1*(1-EMAM1)+daily_prices[MAcounter1-1]*EMAM1
            WMA1 += daily_prices[MAcounter1-1]*MAcounter1
            MAcounter1 += 1
        while(MAcounter2<=MAlength2):
            SMA2 += daily_prices[current_day-MAcounter2]
            EMA2 = EMA2*(1-EMAM2)+daily_prices[MAcounter2-1]*EMAM2
            WMA2 += daily_prices[MAcounter2-1]*MAcounter2
            MAcounter2 += 1
        #this loop runs through all the days in the daily prices list
        while(current_day<len(daily_prices)):
            #these two loops set the moving average values.
            #FOR WMAs: it does some cool things. it adds the properly weighted next day, then subtracts the sum of all the days in the average equally weighted to bring down the weight of everything else by 1. i.e. 1*day1+2*day2+3*day3 becomes 1*day1+2*day2+3*day3+3*day4-(day1+day2+day3) which gets to the new average
            WMA1 += MAlength1*daily_prices[current_day]-SMA1
            WMA1 += MAlength2*daily_prices[current_day]-SMA2
            #FOR SMAs: it adds the next day and subtracts the last day
            SMA1 += (daily_prices[current_day]-daily_prices[current_day-MAlength1])
            SMA2 += (daily_prices[current_day]-daily_prices[current_day-MAlength2])
            #FOR EMAs: weights down the current average and adds the weighted new day
            EMA1 = EMA1*(1-EMAM1)+daily_prices[current_day]*EMAM1
            EMA2 = EMA2*(1-EMAM2)+daily_prices[current_day]*EMAM2
            #this part will do the simulated trading for each strategy, in the order simple vs simple, exponential vs exponential, weighted vs weighted, simple vs exponential, exponential vs weighted, simple vs weighted
            #svs
            if((SMA1/MAlength1)>(SMA2/MAlength2) and cash[0] == True):
                cash[0] = False
                ending_capital[0] = (ending_capital[0]-fee)/daily_prices[current_day]
            elif((SMA1/MAlength1)<(SMA2/MAlength2) and cash[0] == False):
                cash[0] = True
                ending_capital[0] = ending_capital[0]*daily_prices[current_day]-fee
            #eve
            if(EMA1>EMA2 and cash[1] == True):
                cash[1] = False
                ending_capital[1] = (ending_capital[1]-fee)/daily_prices[current_day]
            elif(EMA1<EMA2 and cash[1] == False):
                cash[1] = True
                ending_capital[1] = ending_capital[1]*daily_prices[current_day]-fee
            #wvw
            if((WMA1/(MAlength1*(MAlength1+1)))>(WMA2/(MAlength2*(MAlength2+1))) and cash[2] == True):
                cash[2] = False
                ending_capital[2] = (ending_capital[2]-fee)/daily_prices[current_day]
            elif((WMA1/(MAlength1*(MAlength1+1)))>(WMA2/(MAlength2*(MAlength2+1))) and cash[2] == False):
                cash[2] = True
                ending_capital[2] = ending_capital[2]*daily_prices[current_day]-fee
            #sve
            if((SMA1/MAlength1)>EMA1 and cash[3] == True):
                cash[3] = False
                ending_capital[3] = (ending_capital[3]-fee)/daily_prices[current_day]
            elif((SMA1/MAlength1)<EMA1 and cash[3] == False):
                cash[3] = True
                ending_capital[3] = ending_capital[3]*daily_prices[current_day]-fee
            #evw
            if(EMA1>(WMA1/(MAlength1*(MAlength1+1))) and cash[4] == True):
                cash[4] = False
                ending_capital[4] = (ending_capital[4]-fee)/daily_prices[current_day]
            elif(EMA1<(WMA1/(MAlength1*(MAlength1+1))) and cash[4] == False):
                cash[4] = True
                ending_capital[4] = ending_capital[4]*daily_prices[current_day]-fee
            #svw
            if((SMA1/MAlength1)>(WMA1/(MAlength1*(MAlength1+1))) and cash[5] == True):
                cash[5] = False
                ending_capital[5] = (ending_capital[5]-fee)/daily_prices[current_day]
            elif((SMA1/MAlength1)<(WMA1/(MAlength1*(MAlength1+1))) and cash[5] == False):
                cash[5] = True
                ending_capital[5] = ending_capital[5]*daily_prices[current_day]-fee
            current_day += 1
        #this simulates selling all your stock on the last day if you are currently invested
        if(cash[0] == False):
            cash[0] = True
            ending_capital[0] = ending_capital[0]*daily_prices[current_day-1]-fee
        if(cash[1] == False):
            cash[1] = True
            ending_capital[1] = ending_capital[1]*daily_prices[current_day-1]-fee
        if(cash[2] == False):
            cash[2] = True
            ending_capital[2] = ending_capital[2]*daily_prices[current_day-1]-fee
        if(cash[3] == False):
            cash[3] = True
            ending_capital[3] = ending_capital[3]*daily_prices[current_day-1]-fee
        if(cash[4] == False):
            cash[4] = True
            ending_capital[4] = ending_capital[4]*daily_prices[current_day-1]-fee
        if(cash[5] == False):
            cash[5] = True
            ending_capital[5] = ending_capital[5]*daily_prices[current_day-1]-fee
        #this part adds the result to the end of the list, then changes the length of the moving averages being tested to check aaaaaaallllll the combos
        results.append([ending_capital[0],[MAlength1,MAlength2],"svs"])
        results.append([ending_capital[1],[MAlength1,MAlength2],"eve"])
        results.append([ending_capital[2],[MAlength1,MAlength2],"wvw"])
        results.append([ending_capital[3],[MAlength1,MAlength2],"sve"])
        results.append([ending_capital[4],[MAlength1,MAlength2],"evw"])
        results.append([ending_capital[5],[MAlength1,MAlength2],"svw"])
        ending_capital = [starting_capital,starting_capital,starting_capital,starting_capital,starting_capital,starting_capital]
        if(test_phase == 1):
            MAlength2 += 1
            if(MAlength2>maximum_MA_time_frame):
                print(len(results))
                MAlength1 += 1
                MAlength2 = MAlength1 + 1
            if(MAlength1 == maximum_MA_time_frame):
                test_phase = 2
                print("halftime")
                MAlength1 = 2
                MAlength2 = 1
        elif(test_phase == 2):
            MAlength1 += 1
            if(MAlength1>maximum_MA_time_frame):
                print(len(results))
                MAlength2 += 1
                MAlength1 = MAlength2 + 1
            if(MAlength2 == maximum_MA_time_frame):
                test_phase = 3
    #calculates how much money would have been made buying on the first day and selling on the last day
    default = starting_capital/daily_prices[0]*daily_prices[len(daily_prices)-1]
    #sorts the results by their return in nondecreasing order and returns them, along with the default value at the end
    sorted_results = sorted(results)
    #adds the default entry to the sorted results list so you can always see it when it prints
    sorted_results.append([default,[0,0],"Default"])
    return[sorted_results]
#i thought this one up to help me short tesla in january 2020
#basically you put an array of a stock's daily prices over the course of a (preferably very long, like a year) rally. then what it does is keeps track of the highest price ever, and then when a stock's price drops from the highest ever, it calculates by how much, and holds that percentage value if it's the largest drop so far. at the end it spits out the largest percentage drop(from the maximum price) that the stock has experienced over the data set. the idea is so that you can then set appropriate stop losses or short selling orders when the stock has a LARGER dip from its maximum price
def largest_percentage_drop(prices):
    count = 0
    highest_price = 0
    low_price = 0
    #largest drop is expressed as the ratio of the low price/highest price, and so will be 0.5 for a 50% drop, 0.9 for a 10% drop, etc.
    largest_drop = 1
    while(count<len(prices)):
        if(prices[count]>highest_price):
            highest_price = prices[count]
        elif((prices[count]/highest_price)<largest_drop):
            largest_drop = prices[count]/highest_price
            low_price = prices[count]
        count += 1
    return[largest_drop,low_price]
#quicksorts an array
def quicksort(array):
    counter = 0
    while(len(array) >= 2 and array[len(array)-1] == array[0] and counter < len(array)):
        array.insert(0,array[len(array)-1])
        array.pop()
        counter += 1
    if(len(array) < 2 or array[len(array)-1]==array[0]):
        return[array]
    else:
        arraylow = []
        arrayhigh = []
        #selects pivot point as rightmost element
        pivot = array[len(array)-1]
        #splits up the array into two smaller arrays
        for x in array:
            if(x<=pivot):
                arraylow.append(x)
            else:
                arrayhigh.append(x)
        arraylow.pop()
        arraylow = quicksort(arraylow)[0]
        arrayhigh.insert(0,pivot)
        arrayhigh = quicksort(arrayhigh)[0]
        arraysorted = arraylow + arrayhigh
        return[arraysorted]
#calculates the mean, median, and mode of a dataset, and returns them all in an array format
def mean(array):
    mean = 0
    for x in array:
        mean += x
    mean = mean/len(array)
    return[mean]
def median(array):
    array1 = quicksort(array)[0]
    if(len(array1)%2 == 0):
        median = (array1[len(array1)//2-1]+array1[len(array1)//2-2])/2
    else:
        median = array1[(len(array1)-1)//2]
    return[median]
def mode(array):
    mode = "there is no mode"
    max_occurences = 0
    occurences = 0
    index = 1
    array1 = quicksort(array)[0]
    while(index<len(array1)):
        if(array1[index] == array1[index-1]):
            occurences += 1
        else:
            occurences = 0
        if(occurences > max_occurences):
            max_occurences = occurences
            mode = array1[index]
        elif(occurences == max_occurences):
            mode = "there is no mode"
        index += 1
    return[mode]
#calculates standard deviation and variance for a dataset, sample is a boolean variable, which is true if you want to take sample standard deviation, and false if you want to take population standard deviation
def standard_deviation_variance(array,sample):
    mean1 = mean(array)[0]
    variance = 0
    for x in array:
        variance += pow((x-mean1),2)
    if(sample):
        variance = variance/(len(array)-1)
    else:
        variance = variance/len(array)
    standard_deviation = math.sqrt(variance)
    return[standard_deviation,variance]
def interquartile_range(array):
    array1 = quicksort(array)[0]
    quartile_size = (len(array1)//2+1)/2
    if(quartile_size%1 == 0):
        quartile1 = array1[int(quartile_size-1)]
        quartile3 = array1[int(len(array)-quartile_size)]
    else:
        quartile1 = (array1[int(quartile_size-0.5)]+array1[int(quartile_size-1.5)])/2
        quartile3 = (array1[int(len(array)-quartile_size+0.5)]+array1[int(len(array)-quartile_size-0.5)])/2
    return[quartile3-quartile1]
#sample is a boolean variable, which is true if you want to take sample standard deviation, and false if you want to take population standard deviation
def z_score(array,point,sample):
    mean1 = mean(array)[0]
    standard_deviation = standard_deviation_variance(array,sample)[0]
    zscore = (point-mean1)/standard_deviation
    return[zscore]
#sample is a boolean variable, which is true if you're working with sample data, false if entire population data
def correlation_coefficient(points,sample):
    x_list = []
    y_list = []
    for i in points:
        x_list.append(i[0])
        y_list.append(i[1])
    meanx = mean(x_list)[0]
    meany = mean(y_list)[0]
    covariance = 0
    for j in points:
        covariance += (j[0]-meanx)*(j[1]-meany)
    if(sample):
        covariance = covariance/(len(points)-1)
    else:
        covariance = covariance/len(points)
    standard_deviation_x = standard_deviation_variance(x_list,sample)[0]
    standard_deviation_y = standard_deviation_variance(y_list,sample)[0]
    r = covariance/(standard_deviation_x*standard_deviation_y)
    return[r,covariance]
#put component as true if vectors are in coordinate form, false if it's in polar form, fill in the variables accordingly. fill in z with 0 if you're working in 2D
def dot_product_2D(component,x_or_magnitude1,y_or_degrees1,z1,x_or_magnitude2,y_or_degrees2,z2):
    if(component == True):
        product = x_or_magnitude1*x_or_magnitude2+y_or_degrees1*y_or_degrees2+z1*z2
    else:
        product = x_or_magnitude1*x_or_magnitude2*math.cos(abs(y_or_degrees1-y_or_degrees2)*math.pi/180)
    return product
#vector1 and vector2 must both be 1D lists of real numbers of the SAME LENGTH
def dot_product_nD(vector1, vector2):
    count = 0
    product = 0
    while(count<len(vector1)):
        product += vector1[count]*vector2[count]
        count += 1
    return product
def cross_product(x1,y1,z1,x2,y2,z2):
    x = y1*z2-y2*z1
    y = z1*x2-x1*z2
    z = x1*y2-x2*y1
    return[x,y,z]
#input the coordinates of 2 3d vectors and this will return the coordinates of a unit vector(magnitude of 1) that is perpendicular to the input vectors
def perpendicular_unit_vector(x1,y1,z1,x2,y2,z2):
    #y the ratio of y in terms of z
    constant1 = ((z1*x2)/x1-z2)/(y2-(y1*x2)/x1)
    #this is to sub z into the equation x^2+y^2+z^=1, because it's a unit vector.
    constant2 = ((1+math.pow(constant1,2))+math.pow(((y1*constant1+z1)/x1),2))
    zperp = math.sqrt(1/constant2)
    yperp = constant1*zperp
    xperp = -(y1*yperp+z1*zperp)/x1
    return[xperp,yperp,zperp]
#outputs in degrees
def angle_between_vectors(x1,y1,z1,x2,y2,z2):
    product = dot_product(True,x1,y1,z1,x2,y2,z2)
    magnitude1 = distance([0,0,0],[x1,y1,z1])
    magnitude2 = distance([0,0,0],[x2,y2,z2])
    angle = (180/math.pi)*math.acos(product/(magnitude1*magnitude2))
    return angle
#input two integers, returns their greatest common factor(GCF/GCD)
def greatest_common_divisor(a,b):
    if(b>a):
        c = a
        b = c
        a = b
    c = a%b
    while(c>0):
        a = b
        b = c
        c = a%b
    return b
#this will take 3 of the starting value(a), common ratio(r), number of terms(n), and sum(sum) of a geometric series and find the other one.
#NOTE: put n as "inf" if it's an infinite series and put None for the value you want found
#UNFINISHED, doesn't work for infinite series or unkown r values
def geometric_series(a,r,n,sum):
    if(n == "inf"):
        return[None]
    else:
        if(sum == None):
            sum = a*((1-math.pow(r,n))/(1-r))
        elif(a == None):
            a = sum/((1-math.pow(r,n))/(1-r))
        elif(n == None):
            n = math.log((a+r*sum-sum)/a,r)
        elif(r == None):
            return[None]
    return [a,r,n,sum]
#will approximate a definite integral using simpsons rule
#NOTE: put the operations of the function IN the METHOD in this method
#NOTE 2: n has to be and even number. don't hit me with the n = 75 or its gonna break
#NOTE 3: a<b
def integral_simpsons_rule(n,a,b):
    #insert f(x) in this method, i.e. sin((e^x)^2) or something
    def f(x):
        y = math.pow(math.cos(x),4)*math.sin(x)
        return y
    dx = (b-a)/n
    area = f(a)
    count = 1
    while(count<n):
        area += 4*f(a+count*dx)
        area += 2*f(a+(count+1)*dx)
        count += 2
    area -= f(a+(count-1)*dx)
    area = (area*dx)/3
    return[area]
def nchoosek(n,k):
    return[math.factorial(n)/(math.factorial(k)*math.factorial(n-k))]
#returns the complex quotient (ra+i*(ia))/(rb+i*(ib))
def complex_division(ra, ia, rb, ib):
    denom = rb*rb+ib*ib
    rnum = ra*rb+ia*ib
    inum = ia*rb-ra*ib
    return [rnum/denom,inum/denom]
print(complex_division(6,-8,2,-1))