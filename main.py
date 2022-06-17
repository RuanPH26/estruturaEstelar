import numpy as np


#
#  statstar.py
#
#  This program is a Python version of the code STATSTAR written
#  in FORTRAN and described in Carroll and Ostlie.  Program
#  has been modified and converted to Python in Fall 2010
#  by Prof. D. Schiminovich at Columbia University
#  for use in C3101 Stellar Structure and Evolution
#
#  Notes from FORTRAN code:
#
#  This program is listed in Appendix H of both "An Introduction
#  to Modern Astrophysics," Bradley W. Carroll and Dale A. Ostlie,
#  Addison-Wesley Publishing Company, copyright 1996, and "An
#  Introduction to Modern Stellar Astrophysics," Dale A. Ostlie
#  and Bradley W. Carroll, Addison-Wesley Publishing Company,
#  copyright 1996.
#
#  This program will calculate a static stellar model using the
#  equations developed in the text.  The user is expected to supply the
#  star's mass, luminosity, effective temperature, and composition
#  (X and Z).  If the choices for these quantities are not consistent
#  with the central boundary conditions, an error message will be
#  generated and the user will then need to supply a different set of
#  initial values.

#  SUBROUTINES

# A sub-rotina STARTMDL calcula os valores de M_r, L_r, P e T, próxima a
# superfície da estrela usando expansões da estrutura estelar
# equações (M_r e L_r são considerados constantes).
#
def STARTMDL(i, deltar, X, Z, mu, Rs, Te, r_i, M_ri, L_ri, tog_bf, irc, cst):
    #
    # X = fração de massa do hidrogênio
    # Z = fração de massa do metal
    # mu = peso molecular médio
    # Rs = raio da estrela
    # r_i = raio da casca
    # M_ri = massa dentro da casca neste raio
    # (assumida constante no início do modelo)
    # L_ri = luminosidade fechada, assumida constante
    #
    r = r_i + deltar
    M_rip1 = M_ri
    L_rip1 = L_ri
    #
    #  Esta é a aproximação radiativa (despreze a pressão de radiação
    #  e opacidade de espalhamento de elétrons)
    #  veja Prialnik Eq. 5.1, 5.3 e Sec. 3.7 ou C&O Eqs. (H.1), (H.2), (9.19),
    #  e (9.20).
    #
    if (irc == 0):
        T_ip1 = cst.G * M_rip1 * mu * cst.m_H / (4.25e0 * cst.k_B) * (1.0e0 / r - 1.0e0 / Rs)
        A_bf = 4.34e25 * Z * (1.0e0 + X) / tog_bf
        A_ff = 3.68e22 * cst.g_ff * (1.0e0 - Z) * (1.0e0 + X)
        Afac = (A_bf + A_ff)
        P_ip1 = np.sqrt((1.0e0 / 4.25e0) * (16.0e0 / 3.0e0 * np.pi * cst.a * cst.c) * (cst.G * M_rip1 / L_rip1) * (
                    cst.k_B / (Afac * mu * cst.m_H))) * T_ip1 ** 4.25e0

    #
    #  Esta é a aproximação convectiva# veja Prialnik Sec 6.5, 6.6 ou C&O Eqs. (H.3) e (10,75).
    #
    else:
        T_ip1 = cst.G * M_rip1 * mu * cst.m_H / cst.k_B * (1.0e0 / r - 1.0e0 / Rs) / cst.gamrat
        P_ip1 = cst.kPad * T_ip1 ** cst.gamrat

    return (r, P_ip1, M_rip1, L_rip1, T_ip1)


#
# A subrotina EOS calcula os valores de densidade, opacidade,
# razão do fator guilhotina-gaunt e a taxa de geração de energia em
# um raio r.
#
def EOS(X, Y, Z, XCNO, mu, P, T, izone, cst, idrflg):
    #
    #  Resolva a densidade da lei do gás ideal (remova a pressão de radiação)
    #  veja Prialnik Eq. 5.5 ou C&O Eq. (10.26).
    #
    if ((T < 0.0e0) or (P < 0.0e0)):
        print(' Alguma coisa está um pouco errada aqui.')
        print(' Você está me pedindo para lidar com uma temperatura negativa')
        print(' ou uma pressão negativa. sinto muito mas isso não está no meu')
        print(' contrato! Você terá que tentar novamente com diferentes')
        print(' condições iniciais.')
        print(' Caso ajude, detectei o problema na zona ', izone)
        print(' com as seguintes condições:')
        print('T = ', T, ' K')
        print('P = ', P, ' dynes/cm**2')
        return (0.0, 0.0, 0.0, 0.0, 1)

    Prad = cst.a * T ** 4 / 3.0e0
    Pgas = P - Prad

    rho=(mu*cst.m_H/cst.k_B)*(Pgas/T)

    if (idrflg == 1 and Pgas <= 0):
        Pgas = 5000 * P
    elif (idrflg == 2 and Pgas <= 0):
        Pgas = 0.0001 * P
    elif (idrflg == 0 and Pgas <= 0):
        Pgas = P
    if (T > 0. and Pgas > 0.):
        rho = (mu * cst.m_H / cst.k_B) * (Pgas / T)

    if (rho < 0.0e0):
        print(' Desculpe, mas foi detectada uma densidade negativa.')
        print(' minha rotina de equação de estado está um pouco confusa com esse novo')
        print(' sistema físico que você criou. A pressão de radiação')
        print(' é provavelmente muito grande, o que implica que a estrela é instável.')
        print(' Por favor, tente algo um pouco menos radical da próxima vez.')
        print(' Caso ajude, detectei o problema na zona ', izone)
        print(' com as seguintes condições:')
        print('T       = {0:12.5E} K'.format(T))
        print('P_total = {0:12.5E} dynes/cm**2'.format(P))
        print('P_rad   = {0:12.5E} dynes/cm**2'.format(Prad))
        print('P_gas   = {0:12.5E} dynes/cm**2'.format(Pgas))
        print('rho     = {0:12.5E} g/cm**3'.format(rho))
        return (rho, 0.0, 0.0, 0.0, 1)

    #
    # Calcular a opacidade, incluindo a razão do fator guilhotina para o gaunt#
    # ver Novotny (1973), p. 469. k_bf, k_ff e k_e são os sem limite,
    # livre-livre e opacidades de espalhamento de elétrons, dadas por Prialnik Sec 3.7 ou C&O Eqs. (9.19),
    # (9.20) e (9.21), respectivamente.
    #

    tog_bf = 2.82e0 * (rho * (1.0e0 + X)) ** 0.2e0
    k_bf = 4.34e25 / tog_bf * Z * (1.0e0 + X) * rho / T ** 3.5e0
    k_ff = 3.68e22 * cst.g_ff * (1.0e0 - Z) * (1.0e0 + X) * rho / T ** 3.5e0
    k_e = 0.2e0 * (1.0e0 + X)

    if ((T > 3000. and T < 6000.) and (rho > 1E-10 and rho < 1E-5) and (Z > 0.001 and Z < 0.03)):
        k_Hminus = 2.498e-31 * (Z / 0.02) * np.sqrt(rho) * T ** 9  # Eq. (9.28) C&O in cgs
    else:
        k_Hminus = 0

    kappa = (k_bf + k_ff + k_e + k_Hminus)

    #
    # Calcular a geração de energia pela cadeia pp e o ciclo CNO. Esses
    # são calculados usando as Eqs. (10,49) e (10,53) Prialnik Eq. 4,6, 4,7 + fator de triagem, que vêm de
    # Fowler, Caughlan e Zimmerman (1975). O fator de triagem para o
    # cadeia pp é calculada como fpp# ver Clayton (1968), p. 359ss.
    #
    #     Calculate 1/3 and 2/3 for convenience
    oneo3 = 0.333333333e0
    twoo3 = 0.666666667e0

    T6 = T * 1.0e-6
    T8 = T * 1.0e-8

    # Cadeias PP (ver Hansen e Kawaler, Eq. 6.65, 6.73 e 6.74)
    fx = 0.133e0 * X * np.sqrt((3.0e0 + X) * rho) / T6 ** 1.5e0
    fpp = 1.0e0 + fx * X
    psipp = 1.0e0 + 1.412e8 * (1.0e0 / X - 1.0e0) * np.exp(-49.98 * T6 ** ((-1.0) * oneo3))
    Cpp = 1.0e0 + 0.0123e0 * T6 ** oneo3 + 0.0109e0 * T6 ** twoo3 + 0.000938e0 * T6
    epspp = 2.38e6 * rho * X * X * fpp * psipp * Cpp * T6 ** (-twoo3) * np.exp(-33.80e0 * T6 ** (-oneo3))

    # Ciclo CNO (Kippenhahn e Weigert, Eq. 18.65)
    CCNO = 1.0e0 + 0.0027e0 * T6 ** oneo3 - 0.00778e0 * T6 ** twoo3 - 0.000149e0 * T6
    epsCNO = 8.67e27 * rho * X * XCNO * CCNO * T6 ** (-twoo3) * np.exp(-152.28e0 * T6 ** (-oneo3))

    # Queima de hélio (Kippenhahn e Weigert, Eq. 18.67)
    f3a = 1.  # shielding
    eps_He = 5.09e11 * (rho ** 2) * (Y ** 3) / T8 ** 3 * f3a * np.exp(-44.027 / T8)

    epslon = epspp + epsCNO + eps_He

    return (rho, kappa, epslon, tog_bf, 0)


#
#  Os quatro subprogramas de função a seguir calculam os gradientes de
# # pressão, massa, luminosidade e temperatura em r.

def dPdr(r, M_r, rho, cst):
    return (-cst.G * rho * M_r / (r ** 2))


def dMdr(r, rho, cst):
    return (4.0e0 * np.pi * rho * r ** 2)


def dLdr(r, rho, epslon, cst):
    return (4.0e0 * np.pi * rho * epslon * r ** 2)


def dTdr(r, M_r, L_r, T, rho, kappa, mu, irc, cst):
    if (irc == 0):
        return (-(3.0e0/(16.0e0*np.pi*cst.a*cst.c))*kappa*rho/T**3*L_r/r**2)
        #return -(kappa * rho / T ** 3) * (L_r / 4. / np.pi / r ** 2) * (3 / 4 / cst.a / cst.c)
    #  Este é o gradiente de temperatura convectivo adiabático (Prialnik Eq. 6.29 ou C&O Eq. 10.81).
    else:
        return -(1.0 / cst.gamrat) * (mu * cst.m_H / cst.k_B) * (cst.G * M_r / r ** 2)


#
# Runge-kutta algorithm
#
def RUNGE(f_0, dfdr, r_0, deltar, irc, X, Y, Z, XCNO, mu, izone, cst, idrflg):
    f_temp = np.zeros(4)
    f_i = np.zeros(4)

    dr12 = deltar / 2.0e0
    dr16 = deltar / 6.0e0
    r12 = r_0 + dr12
    r_i = r_0 + deltar
    #
    #  Calcular derivadas intermediárias das quatro estelares fundamentais
    #  equações de estrutura encontradas na Sub-rotina FUNDEQ.
    #

    k1 = deltar * dfdr

    df1, ierr = FUNDEQ(f_0, r_0 + deltar / 2, f_0 + k1 / 2, irc, X, Y, Z, XCNO, mu, izone, cst, idrflg)
    if (ierr != 0):
        return (f_i, ierr)
    k2 = deltar * df1

    df2, ierr = FUNDEQ(f_0, r_0 + deltar / 2, f_0 + k2 / 2, irc, X, Y, Z, XCNO, mu, izone, cst, idrflg)
    if (ierr != 0):
        return (f_i, ierr)
    k3 = deltar * df2

    df3, ierr = FUNDEQ(f_0, r_0 + deltar, f_0 + k3, irc, X, Y, Z, XCNO, mu, izone, cst, idrflg)
    if (ierr != 0):
        return (f_i, ierr)
    k4 = deltar * df3
    #
    #  Calcule as variáveis na próxima camada (i + 1).
    #
    f_i = f_0 + (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)

    return (f_i, 0)


#
#  Esta sub-rotina retorna as derivadas necessárias para RUNGE, o
# # Rotina de integração Runge-Kutta.
#
#      Subroutine FUNDEQ(r, f, irc, X, Z, XCNO, mu, izone,cst)

def FUNDEQ(f0, r, f, irc, X, Y, Z, XCNO, mu, izone, cst, idrflg):
    dfdr = np.zeros(4)
    P0, M_r0, L_r0, T0 = f0
    P, M_r, L_r, T = f

    rho, kappa, epslon, tog_bf, ierr = EOS(X, Y, Z, XCNO, mu, P, T, izone, cst, idrflg)

    dlnPdlnT = ((T0 + T) / (P0 + P)) * ((P0 - P) / (T0 - T))

    if (dlnPdlnT < cst.gamrat):
        irc = 1
    else:
        irc = 0

    dfdr[0] = dPdr(r, M_r, rho, cst)
    dfdr[1] = dMdr(r, rho, cst)
    dfdr[2] = dLdr(r, rho, epslon, cst)
    dfdr[3] = dTdr(r, M_r, L_r, T, rho, kappa, mu, irc, cst)

    return (dfdr, ierr)


def StatStar(Msolar, Lsolar, Te, X, Z, do_plots=True):
    #
    # Programa principal para calcular a estrutura estelar
    #
    # Variáveis, parâmetros de tempo de execução e configurações:
    #
    # deltar = etapa de integração do raio
    # idrflg = define o tamanho do sinalizador
    # = 0 (tamanho do passo de superfície inicial de Rs/1000.)
    # = 1 (tamanho padrão do passo de Rs/100.)
    # = 2 (tamanho do passo principal de Rs/5000.)
    #
    # Nstart = número de passos para os quais as equações iniciais devem ser usadas
    # (assume-se que a zona mais externa é radiativa)
    # Nstop = número máximo de zonas permitidas na estrela
    # Igoof = flag de condição do modelo final
    # = -1 (número de zonas excedido; também o valor inicial)
    # = 0 (bom modelo)
    # = 1 (a densidade do núcleo era extrema)
    # = 2 (a luminosidade do núcleo era extrema)
    # = 3 (temperatura central extrapolada é muito baixa)
    # = 4 (a massa tornou-se negativa antes que o centro fosse alcançado)
    # = 5 (a luminosidade tornou-se negativa antes que o centro fosse alcançado)
    # X, Y, Z = frações de massa de hidrogênio, hélio e metais
    # T0, P0 = temperatura e pressão da superfície (T0 = P0 = 0 é assumido)
    # Ms, Ls, Rs = massa, luminosidade e raio da estrela (unidades cgs)
    #

    nsh = 999  # maximum number of shells

    # initialize r, P, M_r, L_r, T, rho, kappa, epslon, dlPdlT to
    # an array with nsh elements

    r = np.zeros(nsh, float)
    P = np.zeros(nsh, float)
    M_r = np.zeros(nsh, float)
    L_r = np.zeros(nsh, float)
    T = np.zeros(nsh, float)
    rho = np.zeros(nsh, float)
    kappa = np.zeros(nsh, float)
    epslon = np.zeros(nsh, float)
    tog_bf = np.zeros(nsh, float)
    dlPdlT = np.zeros(nsh, float)

    # initialize other variables

    deltar = 0.0
    XCNO = 0.0
    mu = 0.0
    Ms = 0.0
    Ls = 0.0
    Rs = 0.0
    T0 = 0.0
    P0 = 0.0
    Pcore = 0.0
    Tcore = 0.0
    rhocor = 0.0
    epscor = 0.0
    rhomax = 0.0
    Rsolar = 0.0

    Y = 1.0 - (X + Z)

    # tog_bf = constante de opacidade sem limites
    # (proporção de guilhotina para fatores magros)
    tog_bf0 = 0.01

    # initialize variables used for 4 structure equation derivatives

    f_im1 = np.zeros(4, float)
    dfdr = np.zeros(4, float)
    f_i = np.zeros(4, float)

    # Various run time parameters

    Nstart = 1
    Nstop = 999
    Igoof = -1
    ierr = 0
    P0 = 0.0
    T0 = 0.0
    dlPlim = 99.9
    debug = 0

    #
    #  Atribuir valores a constantes (unidades cgs)
    #      # Rsun = raio do Sol
    #      # Msun = massa do Sol
    #      # Lsun = luminosidade do Sol
    #      # sigma = constante de Stefan-Boltzmann
    #      #c = velocidade da luz no vácuo
    #      # a = 4*sigma/c (constante de pressão de radiação)
    #      # G = constante gravitacional universal
    #      # K_B = constante de Boltzmann
    #      # m_H = massa do átomo de hidrogênio
    #      # gama = 5/3 (gama adiabático para um gás monoatômico)
    #      # gama = gama/(gama-1)
    #      # kPad = P/T**(gama/(gama-1)) (constante da lei dos gases adiabáticos)
    #      # g_ff = fator de gaunt de opacidade livre-livre (assumido como unidade)
    #

    Rsun = 6.9599e8
    Msun = 1.989e33
    Lsun = 3.826e33

    # Por conveniência, armazenaremos a maioria das constantes em uma única variável

    class Constants:
        pass

    cst = Constants()

    cst.sigma = 5.67051e-5
    cst.c = 2.99792458e10
    cst.a = 7.56591e-15
    cst.G = 6.67259e-8
    cst.k_B = 1.380658e-16
    cst.m_H = 1.673534e-24
    cst.gamma = 5.0e0 / 3
    cst.g_ff = 1.0e0

    #
    #  Selecione a fração de massa CNO para ser 50% de Z.
    #
    XCNO = Z / 2.0e0
    #
    #  Calcule a massa, a luminosidade e o raio da estrela.
    #      # O raio é calculado a partir da Eq de Prialnik. 1.3 (C&O Eq. 3.17).
    #
    Ms = Msolar * Msun
    Ls = Lsolar * Lsun
    Rs = np.sqrt(Ls / (4 * np.pi * cst.sigma * Te ** 4))
    Rsolar = Rs / Rsun
    #
    #  Comece com um tamanho de passo muito pequeno, pois as condições da superfície variam
    #  rapidamente.
    #
    deltar = -Rs / 1000.
    idrflg = 0
    #
    #  Calcule o peso molecular médio mu assumindo ionização completa
    #      # (ver Prialnik Eq. 3.18, 3.25 ou C&O Eq. 10.21).
    #
    mu = 1.0e0 / (2.0 * X + 0.75 * Y + 0.5 * Z)

    #
    #  Calcule o delimitador entre convecção adiabática e radiação
    #  (ver Prialnik Eq. 6.28 ou C&O Eq. 10.87).
    #
    cst.gamrat = cst.gamma / (cst.gamma - 1.0e0)
    #
    # Inicializa os valores de r, P, M_r, L_r, T, rho, kappa e epslon em
    # a superfície. A zona mais externa é considerada a zona 1. A zona
    # número aumenta em direção ao centro.
    #

    initsh = 0
    r[initsh] = Rs
    M_r[initsh] = Ms
    L_r[initsh] = Ls
    T[initsh] = T0
    P[initsh] = P0
    tog_bf[initsh] = tog_bf0

    if (P0 <= 0.0) or (T0 <= 0.0):
        rho[initsh] = 0.0
        kappa[initsh] = 0.0
        epslon[initsh] = 0.0
        tog_bf[initsh] = 0.01
    else:
        rho[initsh], kappa[initsh], epslon[initsh], tog_bf[initsh], ierr = EOS(X, Y, Z, XCNO, mu, P[initsh], T[initsh],
                                                                               0, cst, idrflg)
        if ierr != 0:
            print("estamos parando agora")
            istop = 0

    #
    #  Aplique soluções de superfície aproximadas para iniciar a integração,
    #  supondo transporte de radiação na zona mais externa (o loop do 20).
    #  irc = 0 para radiação, irc = 1 para convecção.
    #  Assumir valores iniciais arbitrários para kPad e dlPdlT.
    #  dlPdlT = dlnP/dlnT (ver Prialnik Eq. 6.28 ou C&O Eq. 10.87)
    #
    cst.kPad = 0.3e0
    irc = 0
    dlPdlT[initsh] = 4.5

    for i in range(0, Nstart):
        ip1 = i + 1
        r[ip1], P[ip1], M_r[ip1], L_r[ip1], T[ip1] = STARTMDL(i, deltar, X, Z, mu, Rs, Te, r[i], M_r[i], L_r[i],
                                                              tog_bf[i], irc, cst)
        rho[ip1], kappa[ip1], epslon[ip1], tog_bf[ip1], ierr = EOS(X, Y, Z, XCNO, mu, P[ip1], T[ip1], ip1, cst, idrflg)

        dlnPdlnT0 = ((T[i] + T[ip1]) / (P[i] + P[ip1])) * ((P[i] - P[ip1]) / (T[i] - T[ip1]))

        if ierr != 0:
            print('Values from the previous zone are:')
            print('r/Rs      = {0:12.5E}'.format(r[i] / Rs))
            print('rho       = {0:12.5E}  g/cm**3'.format(rho[i]))
            print('M_r/Ms    = {0:12.5E}'.format(M_r[i] / Ms))
            print('kappa     = {0:12.5E}  cm**2/g'.format(kappa[i]))
            print('T         = {0:12.5E}  K'.format(T[i]))
            print('epsilon   = {0:12.5E}  ergs/g/s'.format(epslon[i]))
            print('P         = {0:12.5E}  dynes/cm**2'.format(P[i]))
            print('L_r/Ls    = {0:12.5E}'.format(L_r[i] / Ls))
            break
        #
        #  Determine se a convecção estará operando na próxima zona
        #  calculando dlnP/dlnT numericamente entre as zonas i e i+1 [ip1].
        #  Atualize a constante de gás adiabático, se necessário.
        #
        if (i > initsh):
            dlPdlT[ip1] = ((T[i] + T[ip1]) / (P[i] + P[ip1])) * ((P[i] - P[ip1]) / (T[i] - T[ip1]))
        else:
            dlPdlT[ip1] = ((T[i] + T[ip1]) / (P[i] + P[ip1])) * ((P[i] - P[ip1]) / (T[i] - T[ip1]))

        if (dlPdlT[ip1] < cst.gamrat):
            irc = 1
        else:
            irc = 0
            cst.kPad = P[ip1] / T[ip1] ** cst.gamrat

        #
        #  Teste para ver se a suposição de superfície de massa constante ainda é
        #  válida.
        #
        deltaM = deltar * dMdr(r[ip1], rho[ip1], cst)
        M_r[ip1] = M_r[i] + deltaM
        if (np.abs(deltaM) > (1.0e-8 * Ms)):
            if (ip1 > 1):
                ip1 = ip1 - 1
                print(' A variação na massa tornou-se maior que 0,001*Ms')
                print(' deixando o loop de aproximação antes que Nstart fosse alcançado.')
                break
    #
    #
    #  =============  MAIN LOOP ==========================
    #
    #  Este é o loop de integração principal. As suposições de constante
    #  massa e luminosidade interior não são mais aplicadas.
    #
    #  Inicialize a rotina Runge-Kutta especificando as quantidades da zona i-1
    #  e seus derivados. Observe que a pressão, massa, luminosidade,
    #  temperatura são armazenados nos locais de memória f_im1(1),
    #  f_im1(2), f_im1(3) e f_im1(4), respectivamente. Os derivados de
    #  essas quantidades em relação ao raio são armazenadas em dfdr(1),
    #  dfdr(2), dfdr(3) e dfdr(4). Por fim, os valores resultantes para
    #  P, M_r, L_r e T são retornados da rotina Runge-Kutta em
    #  f_i(1), f_i(2), f_i(3) e f_i(4), respectivamente.
    #  # As equações da estrutura estelar
    #  # dPdr (Prialnik 5.1 ou C&O Eq. 10.7)
    #  # dMdr (Prialnik 5.2 ou C&O Eq. 10.8)
    #  # dLdr (Prialnik 5.4 ou C&O Eq. 10.45) e
    #  # dTdr (Prialnik 5.3 ou C&O Eq. 10.61 ou Eq. 10.81)
    #  # são calculados em chamadas de função, definidas anteriormente no código.
    #

    Nsrtp1 = ip1 + 1

    if (ierr != 0):  # sair se chegamos a este ponto após um erro na inicialização
        Nstop = Nsrtp1 - 1
        istop = Nstop

    for i in range(Nsrtp1, Nstop):

        im1 = i - 1
        f_im1[0] = P[im1]
        f_im1[1] = M_r[im1]
        f_im1[2] = L_r[im1]
        f_im1[3] = T[im1]
        dfdr[0] = dPdr(r[im1], M_r[im1], rho[im1], cst)
        dfdr[1] = dMdr(r[im1], rho[im1], cst)
        dfdr[2] = dLdr(r[im1], rho[im1], epslon[im1], cst)
        dfdr[3] = dTdr(r[im1], M_r[im1], L_r[im1], T[im1], rho[im1], kappa[im1], mu, irc, cst)

        f_i, ierr = RUNGE(f_im1, dfdr, r[im1], deltar, irc, X, Y, Z, XCNO, mu, i, cst, idrflg)

        if (ierr != 0):
            print(' O problema ocorreu na rotina Runge-Kutta')
            print(' Os valores da zona anterior são:')
            print('r/Rs    = {0:12.5e}'.format(r[im1] / Rs))
            print('rho     = {0:12.5e} g/cm**3'.format(rho[im1]))
            print('M_r/Ms  = {0:12.5e}'.format(M_r[im1] / Ms))
            print('kappa   = {0:12.5e} cm**2/g'.format(kappa[im1]))
            print('T       = {0:12.5e} K'.format(T[im1]))
            print('epsilon = {0:12.5e} ergs/g/s'.format(epslon[im1]))
            print('P       = {0:12.5e} dynes/cm**2'.format(P[im1]))
            print('L_r/Ls  = {0:12.5e}'.format(L_r[im1] / Ls))
            break
        #
        #  Atualize os parâmetros estelares para a próxima zona, incluindo a adição
        #  dr para o raio antigo (note que dr < 0 já que a integração é
        #  para dentro).
        #
        r[i] = r[im1] + deltar
        P[i] = f_i[0]
        M_r[i] = f_i[1]
        L_r[i] = f_i[2]
        T[i] = f_i[3]
        #
        #  Calcule a densidade, opacidade e taxa de geração de energia para
        #  esta zona.
        #
        rho[i], kappa[i], epslon[i], tog_bf[i], ierr = EOS(X, Y, Z, XCNO, mu, P[i], T[i], i, cst, idrflg)

        if (ierr != 0):
            print(' Os valores da zona anterior são:')
            print('r/Rs    = {0:12.5e}'.format(r[im1] / Rs))
            print('rho     = {0:12.5e} g/cm**3'.format(rho[im1]))
            print('M_r/Ms  = {0:12.5e}'.format(M_r[im1] / Ms))
            print('kappa   = {0:12.5e} cm**2/g'.format(kappa[im1]))
            print('T       = {0:12.5e} K'.format(T[im1]))
            print('epsilon = {0:12.5e} ergs/g/s'.format(epslon[im1]))
            print('P       = {0:12.5e} dynes/cm**2'.format(P[im1]))
            print('L_r/Ls  = {0:12.5e}'.format(L_r[im1] / Ls))
            istop = i
            break

        if (debug == 1): print(i, r[i], M_r[i], L_r[i], T[i], P[i], rho[i], kappa[i], epslon[i], tog_bf[i])

        #
        #  Determine se a convecção estará operando na próxima zona por
        #  calculando dlnP/dlnT e comparando-o com gamma/(gamma-1)
        #  (ver Prialnik Eq. 6.28 ou C&O Eq. 10.87). Defina o sinalizador de convecção adequadamente.
        #
        dlPdlT[i] = np.log(P[i] / P[im1]) / np.log(T[i] / T[im1])
        #dlPdlT[i] = ((T[im1] + T[i])/(P[im1] + P[i]))*((P[im1] - P[i])/(T[im1] - T[i]))

        if (dlPdlT[i] < cst.gamrat):
            irc = 1
        else:
            irc = 0

        #
        #  Verifique se o centro foi alcançado. Se sim, defina Igoof e
        #  estimar as condições centrais rhocor, epscor, Pcore e Tcore.
        #  A densidade central é estimada como sendo a densidade média do
        #  bola central restante, a pressão central é determinada por
        #  usando a expansão de Taylor no centro (Prialnik - Exercício. 5.1; CO Eq. H.4)
        #  e o valor central para a energia
        #  taxa de geração é calculada para ser o interior restante
        #  luminosidade dividida pela massa da bola central. finalmente, o
        #  temperatura central é calculada aplicando a lei dos gases ideais
        #  (onde a pressão de radiação é desprezada).
        #

        ####################################################################
        # use python para extrapolar para o núcleo r[i], Qm, L_r[i], T[i], P[i], rho[i], kappa[i],epslon[i], clim, rcf, dlPdlT[i]
        from scipy import interpolate

        f_Lr = interpolate.interp1d(r[0:i], L_r[0:i], fill_value='extrapolate')
        f_T = interpolate.interp1d(r[0:i], T[0:i], fill_value='extrapolate')
        f_P = interpolate.interp1d(r[0:i], P[0:i], fill_value='extrapolate')
        f_rho = interpolate.interp1d(r[0:i], rho[0:i], fill_value='extrapolate')
        f_eps = interpolate.interp1d(r[0:i], epslon[0:i], fill_value='extrapolate')

        if ((r[i] <= np.abs(deltar)) and ((L_r[i] >= (0.1e0 * Ls)) or (M_r[i] >= (0.01e0 * Ms)))):
            #   Acerte o centro antes que a massa/luminosidade se esgote
            Igoof = 6
        elif (L_r[i] <= 0.0e0):
            #   Obteve luminosidade central negativa
            Igoof = 5
            # rho: Taylor expansion at center
            rhocor = M_r[i] / (4.0e0 / 3.0e0 * np.pi * r[i] ** 3)
            if (M_r[i] != 0.0e0):
                # taxa de geração de energia
                epscor = L_r[i] / M_r[i]
            else:
                epscor = 0.0e0
            # P:  Taylor expansion at center
            Pcore = P[i] + 2.0e0 / 3.0e0 * np.pi * cst.G * rhocor ** 2 * r[i] ** 2
            # Assume ideal gas
            Tcore = Pcore * mu * cst.m_H / (rhocor * cst.k_B)
        elif (M_r[i] <= 0.0e0):
            Igoof = 4  # O modelo tem um buraco no centro (densidade negativa!)
            Rhocor = 0.0e0
            epscor = 0.0e0
            Pcore = 0.0e0
            Tcore = 0.0e0
        elif ((r[i] < (0.02e0 * Rs)) and ((M_r[i] < (0.01e0 * Ms)) and ((L_r[i] < 0.1e0 * Ls)))):
            # se atingimos <2% do raio da estrela,
            # <1% de massa fechada e <10% de luminosidade então....
            # rho: Expansão de Taylor no centro
            #rhocor = M_r[i]/(4./3.*np.pi*r[i]**3)
            # defina a massa de núcleo máxima razoável
            #rhomax = 10.0e0 * (rho[i] / rho[im1]) * rho[i]
            #epscor = L_r[i]/M_r[i]
            ## P: Taylor expansion at center
            #Pcore  = P[i] + 2.0e0/3.0e0*np.pi*cst.G*rhocor**2*r[i]**2
            #   # Assume ideal gas
            #Tcore  = Pcore*mu*cst.m_H/(rhocor*cst.k_B)
            # Em geral, todos estes devem produzir valores
            # que sobem em direção ao centro (mas não muito alto)

            # core values from python extrapolation
            rhocor = f_rho(r[i])
            epscor = f_eps(r[i])
            Pcore = f_P(r[i])
            Tcore = f_T(r[i])

            if ((np.abs(rhocor - rho[i]) / rho[i] > 0.01) or (rhocor > rhomax)):
                # rho está desligado grande ou pequeno
                Igoof = 1
            elif (np.abs(epscor - epslon[i]) / epslon[i] > 0.01):
                # taxa de geração de energia um pouco fora (baixa)
                Igoof = 2
            elif (Tcore < T[i]):
                # Temperatura um pouco fora (baixa)
                Igoof = 3
            else:
                # o número de shells permitidos foi excedido
                Igoof = 0

        if (Igoof != -1):
            istop = i
            break
        #
        #  É hora de mudar o tamanho do passo?
        #
        if ((idrflg == 0) and (M_r[i] < (0.99e0 * Ms))):
            deltar = (-1.0) * Rs / 100
            idrflg = 1

        if ((idrflg == 1) and (deltar >= (0.5 * r[i]))):
            deltar = (-1.0) * Rs / 5000.0e0

            idrflg = 2

        istop = i

    #
    #  Gere mensagens de aviso para as condições centrais.
    #
    #rhocor = M_r[istop]/(4.0e0/3.0e0*np.pi*r[istop]**3)
    #epscor = L_r[istop]/M_r[istop]
    #Pcore  = P[istop] + 2.0e0/3.0e0*np.pi*cst.G*rhocor**2*r[istop]**2
    #Tcore  = Pcore*mu*cst.m_H/(rhocor*cst.k_B)

    rhocor = f_rho(r[istop])
    epscor = f_eps(r[istop])
    Pcore = f_P(r[istop])
    Tcore = f_T(r[istop])

    rhocor = f_rho(0.)
    epscor = f_eps(0.)
    Pcore = f_P(0.)
    Tcore = f_T(0.)

    if (Igoof != 0):
        if (Igoof == -1):
            print('Desculpe ser o portador de más notícias, mas...')
            print('Seu modelo tem alguns problemas')
            print('O número de shells permitidos foi excedido')

        if (Igoof == 1):
            print('Parece que você está chegando perto,')
            print('no entanto, ainda existem alguns pequenos erros')
            print('A densidade do núcleo parece um pouco fora,')
            print(' densidade deve aumentar suavemente em direção ao centro.')
            print(' A densidade da última zona calculada foi rho = ', rho[istop], ' gm/cm**3')
            print(rhocor, rhomax)
        if (rhocor > 1e10):
            print('Parece que você vai precisar de um degenerado')
            print(' gás de nêutrons e relatividade geral')
            print(' para resolver este núcleo. Quem você pensa que eu sou, Einstein?')

        if (Igoof == 2):
            print('Parece que você está chegando perto,')
            print('no entanto, ainda existem alguns pequenos erros')
            print('O núcleo épsilon parece um pouco fora,')
            print(' epsilon deve variar suavemente perto do centro.')
            print(' O valor calculado para a última zona foi eps =', epslon[istop], ' ergs/g/s')
            print(epscor)
        if (Igoof == 3):
            print('Parece que você está chegando perto,')
            print('no entanto, ainda existem alguns pequenos erros')
            print('Sua temperatura central extrapolada é muito baixa')
            print(' um pouco mais de ajuste fino deve fazê-lo.')
            print(' O valor calculado para a última zona foi T = ', T[istop], ' K')

        if (Igoof == 4):
            print('Desculpe ser o portador de más notícias, mas...')
            print('Seu modelo tem alguns problemas')
            print('Você criou uma estrela com um buraco no centro!')

        if (Igoof == 5):
            print('Desculpe ser o portador de más notícias, mas...')
            print('Seu modelo tem alguns problemas')
            print('Esta estrela tem uma luminosidade central negativa!')

        if (Igoof == 6):
            print('Desculpe ser o portador de más notícias, mas...')
            print('Seu modelo tem alguns problemas')
            print('Você atinge o centro antes da massa e/ou ')
            print('a luminosidade se esgotar!')
    else:
        print('PARABÉNS, ACHO QUE VOCÊ ENCONTROU!')
        print('No entanto, certifique-se de olhar para o seu modelo com cuidado.')

    #
    #  Imprima as condições centrais. Se necessário, estabeleça limites para a
    #  raio central, massa e luminosidade se necessário, para evitar formato
    #  campo estourou.
    #

    Rcrat = r[istop] / Rs
    if (Rcrat < -9.999e0): Rcrat = -9.999e0
    Mcrat = M_r[istop] / Ms
    if (Mcrat < -9.999e0): Mcrat = -9.999e0
    Lcrat = L_r[istop] / Ls
    if (Lcrat < -9.999e0): Lcrat = -9.999e0

    f = open('starmodl_py.dat', 'w')

    f.write('Um Modelo Homogêneo de Sequência Principal\n')
    f.write(' As condições da superfície são:       As condições centrais são:\n')
    f.write(' Mtot = {0:13.6E} Msun          Mc/Mtot     = {1:12.5E}\n'.format(Msolar, Mcrat))
    f.write(' Rtot = {0:13.6E} Rsun          Rc/Rtot     = {1:12.5E}\n'.format(Rsolar, Rcrat))
    f.write(' Ltot = {0:13.6E} Lsun          Lc/Ltot     = {1:12.5E}\n'.format(Lsolar, Lcrat))
    f.write(' Teff = {0:13.6E} K             Density     = {1:12.5E}\n'.format(Te, rhocor))
    f.write(' X    = {0:13.6E}               Temperature = {1:12.5E}\n'.format(X, Tcore))
    f.write(' Y    = {0:13.6E}               Pressure    = {1:12.5E} dynes/cm**2\n'.format(Y, Pcore))
    f.write(' Z    = {0:13.6E}               epsilon     = {1:12.5E} ergs/s/g\n'.format(Z, epscor))
    f.write('                                    dlnP/dlnT   = {0:12.5E}\n'.format(dlPdlT[istop]))

    f.write('Notes:\n')
    f.write(' (1) A massa é listada como Qm = 1.0 - M_r/Mtot, quando Mtot = {0:13.6}\n'.format(Msun))
    f.write(' (2) As zonas convectivas são indicadas por c, as zonas radiativas por r\n')
    f.write(' (3) dlnP/dlnT pode ser limitado a +99.9 ou -99.9# if so it is\n')
    f.write(' labeled by *\n')

    if do_plots:
        import matplotlib.pyplot as plt

        plt.plot(r[0:istop + 1], L_r[0:istop + 1] / L_r[0:istop + 1].max(), label='luminosidade')
        plt.plot(r[0:istop + 1], T[0:istop + 1] / T[0:istop + 1].max(), label='temperatura')
        plt.plot(r[0:istop + 1], M_r[0:istop + 1] / M_r[0:istop + 1].max(), label='massa')
        plt.xlabel('Raio')
        plt.legend()
        plt.show()

    #
    # Imprima dados do centro da estrela para fora, rotulando convectivo
    # ou zonas radiativas por c ou r, respectivamente. Se abs(dlnP/dlnT)
    # excede 99,9, defina um sinalizador de aviso de impressão (*) e defina o limite de saída
    # para +99,9 ou -99,9 conforme apropriado para evitar estouros de campo de formato.
    #
    f.write('   r        M_r       L_r       T        P        rho      kap      eps     dlPdlT\n')

    for ic in range(0, istop + 1):
        i = ic  # istop+1 - ic
        Qm = 1.0e0 - M_r[i] / Ms  # Fração de massa total até o raio

        if (dlPdlT[i] < cst.gamrat):
            rcf = 'c'
        else:
            rcf = 'r'
        if (np.abs(dlPdlT[i]) > dlPlim):
            dlPdlT[i] = np.copysign(dlPlim, dlPdlT[i])
            clim = '*'
        else:
            clim = ' '
        s = '{0:7.3E} {1:7.3E} {2:7.3E} {3:7.3E} {4:7.3E} {5:7.3E} {6:7.3E} {7:6.3E}{8:1s}{9:1s} {10:5.1f}\n'.format(
            r[i], M_r[i], L_r[i], T[i], P[i], rho[i], kappa[i], epslon[i], clim, rcf, dlPdlT[i])
        f.write(s)

    #     Output to screen
    print
    print('***** A integração foi concluída *****')
    print('      O modelo foi armazenado em starmodl_py.dat')
    print

    resdidual = np.mean(
        [np.abs(rhocor - rho[i]) / rho[i], np.abs(epscor - epslon[i]) / epslon[i], np.abs(Tcore - T[i]) / T[i]])
    return Igoof, ierr, istop, resdidual

def main():

# #
# #  Enter desired stellar parameters
# #


       getinp=1  # read in input
       if (getinp == 1):
            Msolar=float(input(' Digite a massa da estrela (em unidades solares):'))
            Lsolar=float(input(' Digite a luminosidade da estrela (em unidades solares):'))
            Te=float(input(' Insira a temperatura efetiva da estrela (em K):'))
            Y=-1.0
            while (Y < 0.0):
                X=float(input(' Digite a fração de massa de hidrogênio (X):'))
                Z=float(input(' Digite a fração de massa de metais (Z):'))
                Y = 1.e0 - X - Z
                if Y < 0:
                    print('Você deve ter X + Z <= 1. Por favor, reinsira a composição.')

       Igoof,ierr,istop, resdidual=StatStar(Msolar,Lsolar,Te,X,Z)
main()

