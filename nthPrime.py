import math
import sys
import random
import time

sys.setrecursionlimit(10000)


counts = {
    'lehmer': 0,
    'primality': 0
}
primeCache = []
smallPiLimit = 1000000
lastSieveSize = smallPiLimit + 1



def estimateNthPrime(n):

    guess = math.log(n) + math.log(math.log(n)) - 1
    guess += (math.log(math.log(n)) - 2) / math.log(n)
    guess -= (pow((math.log(math.log(n))), 2) - (6 * math.log(math.log(n))) + 11) / (2*pow(math.log(n), 2))
    guess += 1 / (pow(math.log(n), 2))
    guess *= n
    guess = round(guess)
    if (guess % 2 == 0):
        guess += 1

    return guess



def smallPrime(n):
    global primeCache
    global lastSieveSize
    while (len(primeCache) <= n):
        primeCache = primeSieve(lastSieveSize)
        lastSieveSize = lastSieveSize * 2

    return primeCache[n]



def primeSieve(n):
    n = n + 1
    # Eratosthenes algorithm to find all primes under n
    array = []
    upperLimit = int(math.ceil(math.sqrt(n)))
    output = [1]

    # Make an array from 2 to (n - 1)
    for i in xrange(0, n):
        array.append(True)

    # Remove multiples of primes starting from 2, 3, 5,...
    for i in xrange(2, upperLimit):
        if (array[i]):
            for j in xrange(i*i, n, i):
                array[j] = False

    # All array[i] set to True are primes
    for i in xrange(2, n):
        if(array[i]):
            output.append(i)

    return output




phiCache = {}

def phi(x, a):

    # print 'phi: ' + str(x) + ' ' + str(a)
    a = int(a)
    x = int(x)
    key = (x, a)
    # key = String(x) + ',' + String(a)

    if (key in phiCache):
        return phiCache[key]

    if (a < 1.1 and a > 0.9):
        return int(math.floor((x+1)/2))

    result = phi(x, a-1) - phi(math.floor(x / smallPrime(a)), a-1)
    phiCache[key] = result
    return result


piCache = {}

def smallPi(x):
    global lastSieveSize
    global primeCache
    global piCache

    if x in piCache:
        return piCache[x]

    if (len(primeCache) == 0):
        primeCache = primeSieve(lastSieveSize)
        lastSieveSize = lastSieveSize * 2


    largest = primeCache[len(primeCache) - 1]
    while (x > largest):
        primeCache = primeSieve(lastSieveSize)
        lastSieveSize = lastSieveSize * 2
        largest = primeCache[len(primeCache) - 1]

    l = 0
    r = len(primeCache) - 1
    while (l != r):
        m = int(math.floor((l + r) / 2))
        if (primeCache[m] > x):
            r = m
        elif (primeCache[m] <= x):
            l = m + 1

    if (primeCache[m] <= x):
        m += 1

    piCache[x] = m - 1
    return m - 1


# TODO: 1 too little for 10^6
def lehmerPi(x):

    # // TODO: Cache bottom N?
    if (x < smallPiLimit):
        return smallPi(x)

    a = lehmerPi(int(math.floor(pow(x, float(1)/float(4)))))
    b = lehmerPi(int(math.floor(pow(x, float(1)/float(2)))))
    c = lehmerPi(int(math.floor(pow(x, float(1)/float(3)))))

    sum = phi(x,a) + (((b+a-2) * (b-a+1)) / 2)
    for i in xrange(a+1, b+1):

        w = int(math.floor(x/smallPrime(i)))
        bi = lehmerPi(int(math.floor(math.pow(w, float(1)/float(2)))))
        sum = sum - lehmerPi(w)
        if (i <= c):
            for j in xrange(i, bi+1):
                sum = sum - (lehmerPi(int(math.floor(x/(smallPrime(j)*smallPrime(i)) ))) - (j - 1))

    return sum


def primeTest(n, k):
    counts['primality'] += 1
    if (n == 2 or n == 3):
        return True
    if (n % 2 == 0 or n < 2):
        return False

    # // Write (n - 1) as 2^s * d
    s = 0
    d = n - 1
    while (d % 2 == 0):
        d = d / 2
        s += 1

    for a in random.sample(xrange(2, n - 2), k):
        flag = True
        x = pow(a, d, n)
        if (x == 1 or x == (n - 1)):
            continue

        for r in xrange(0, s-1):
            x = pow(x, 2, n)

            if (x == 1):
                return False

            if (x == (n - 1)):
                flag = False
                break

        if (flag):
            return False

    return True



def getNthPrime(n):
    beforeLehmer = time.time()
    guess = estimateNthPrime(n)
    guessPi = lehmerPi(guess)

    nBasedTolerance = math.log(n) * 1000
    tolerance = max(3, nBasedTolerance)
    print 'Using tolerance: ' + str(tolerance)
    firstGuess = guess
    piGap = abs(n - guessPi)
    diff = int(round(piGap*math.log(guess)*2))
    if (guessPi < n):
        l = guess
        r = guess + diff
        lPi = guessPi
        rPi = lehmerPi(r)
    elif (guessPi > n):
        r = guess
        l = guess - diff
        lPi = lehmerPi(l)
        rPi = guessPi

    counts['lehmer'] = 2

    while(abs(lPi - rPi) > tolerance):
        m = int(round((l + r)/2))
        mPi = lehmerPi(m)
        counts['lehmer'] += 1

        # // mPi overestimate
        if (mPi >= n):
            rPi = mPi
            r = m
        elif (mPi < n):
            lPi = mPi
            l = m


    guessPi = lPi
    guess = l
    if (guess % 2 == 0):
        guess = guess - 1
        guessPi = lehmerPi(guess)
        counts['lehmer'] += 1

    guess = int(guess)
    print 'Time to do lehmerPi: %s' % (time.time() - beforeLehmer)
    beforePrimeTests = time.time()

    while (guessPi < n):
        guess += 2
        isPrime = primeTest(guess, 30)

        if (isPrime):
            guessPi += 1

    print 'Time to do primality: %s' % (time.time() - beforePrimeTests)
    return guess







start = time.time()
n = int(pow(10, 8)*2);
print 'n: ' + str(n)
print getNthPrime(n)
print 'Execution in:'
print str(time.time() - start)
print counts




