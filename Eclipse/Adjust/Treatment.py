import scipy
import numpy

class Tratamento :
    
    def __init__(self,modelo):

        '''
        Funcao para extrair os transitos individualmente da curva de luz
        
        parâmetro time :: tempo da curva de luz total
        parâmetro flux :: fluxo da curva de luz total
        parâmetro flux_err :: erro do fluxo
        parâmetro u1 :: coeficiente de escurecimento de limbo 1
        parâmetro u2 :: coeficiente de escurecimento de limbo 2
        parâmetro porb :: periodo da órbita em dias 
        parâmetro AU :: semieixo orbital em UA
        parâmetro raioPlan :: raio do planeta em relaçao ao raio da estrela
        parâmetro inc :: angulo de inclinacao em graus 
        parâmetro x0 :: 
        parâmetro nt :: 
        '''
        self.modelo = modelo
        self.u1,self.u2,self.porb,self.time,self.flux,self.flux_err,self.raioPlan,self.AU,self.inc,self.x0,self.nt,self.ts_model, self.mass = modelo.retornaParametros()

    def cut_transit_single(self):
        
        '''
        returns 
        
        parâmetro dur ::  duracao do transito em horas
        parâmetro t_split  :: tempo em horas (igual para todos os transitos)
        parâmetro n_f_split ::curva de luz do transito normalizada
        parâmetro n_f_err_split :: erro do fluxo normalizado
        '''

        if self.u1 == 999:
            self.u1 = 0.59
        if self.u2 == 999:
            self.u2 = 0.0
        
        self.wl, self.u1, self.u2 = self.modelo.limb(0)    
        lc0, ts0 = self.modelo.eclipse_model()
        
        #duração do transito em horas
        
        x = ts0
        y = 1 - lc0
        #yh = 0.002
        yh=max(y)/2.
        kk = numpy.where(y >= yh)
        x1 = min(x[kk])
        x2 = max(x[kk])
        meio = (x1 + x2)/ 2 # valor central dos transitos em fase
        self.dur = 2 * numpy.abs(ts0[min(numpy.where(lc0 < 1))[0]]) #duração total do transito
        # em horas
        
        # dado o valor central de cada transito, determino pontos (+-) a uma distancia de 45%
        # do período. 
        mm = []
        mp = []
        ttt = []

        for i in range(0, int(self.nt)):
            mm.append(self.x0 + self.porb*i - (self.porb*0.45))
            mp.append(self.x0 + self.porb*i + (self.porb*0.45))
            
        # indice que cada valor representa no array do tempo
        for i in range(0, int(self.nt)):
            ttt.append(numpy.where((self.time >= mm[i]) & (self.time <= mp[i])))
        
        # separo os arrays de tempo e fluxo com esses valores. 
        
        self.t_split = []
        f_split = []
        f_err_split = []
        for i in range(0, int(self.nt)):
            self.t_split.append((self.time[ttt[i]]-self.x0-(self.porb*i))*24)   # tempo em horas com o 
            # meio do transito central em zero
            f_split.append(self.flux[ttt[i]])         
            f_err_split.append(self.flux_err[ttt[i]])  
            
        # eliminação do slope
        
        size = [] #quantidade de pontos presentes em até 3* a duração do transito
        for u in range(int(self.nt)):
            size.append(len(numpy.where(numpy.sort(numpy.abs(self.t_split[i])) < self.dur*3.)[0]))  
            
        # considerado dados que apresentam fluxo maior que 0 
        # e dados de tempo com quantidade de pontos com 0.9*numpy.mean(size)
        self.n_f_split = []
        n_f_err_split = []
        for i in range(0, int(self.nt)):
            if len((f_split[i] > 0) & (len(numpy.where(numpy.sort(numpy.abs(self.t_split[i])) < self.dur*3.)[0]) > numpy.mean(size)*.9)):
                ss = numpy.polyfit(self.t_split[i], f_split[i], 1)
                zz = numpy.polyval(ss, self.t_split[i])
                self.n_f_split.append(numpy.array(f_split[i] - zz + 1))
                n_f_err_split.append(f_err_split[i])
            else:
                self.n_f_split.append(f_split[i])
                n_f_err_split.append(f_err_split[i])
                
        #renormalização
        w_flux = []
        for i in range(0, int(self.nt)):
            if len((self.n_f_split[i] > 0) & (len(numpy.where(numpy.sort(numpy.abs(self.t_split[i])) < self.dur*3.)[0]) > numpy.mean(size)*.9)):
                w_flux =  numpy.append(w_flux, self.n_f_split[i])
        m0 = numpy.median(w_flux) 
        for i in range(0, int(self.nt)):
            if len((self.n_f_split[i] > 0) & (len(numpy.where(numpy.sort(numpy.abs(self.t_split[i])) < self.dur*3.)[0]) > numpy.mean(size)*.9)):
                self.n_f_split[i] = self.n_f_split[i]/m0
                n_f_err_split[i] = n_f_err_split[i]/m0

        return self.dur, self.t_split, self.n_f_split, n_f_err_split
        #t_split = tim
        #n_f_split = lcurve
    #--------------------------------------------------#

    def transit_smooth(self,ntransit, selection):
    
        '''
        Funcao para uma curva smooth com n transitos dado um range (selection, onde seleciona os transitos de forma randômica ou os trânsitos mais fundos)
        
        parâmetro ntransit :: numero de transitos para usar na curva smoothed
        parâmetro selection :: 0, usa uma escolha randomica de todos os transitos
        parâmetro se selection :: 1, usa os transitos mais fundos  
        
        returns
        parâmetro time_phased[bb] :: tempo 
        parâmetro smoothed_LC[bb] :: curva de luz Smoothed
        '''
        
        if selection == 0:
            tran_selec = numpy.random.randint(int(self.nt), size=(1, ntransit))[0]
        
        else:
            deepest_transit = []
            for i in range(0, int(self.nt)):
                if len(self.n_f_split[i]) > 0:
                    deepest_transit.append(numpy.mean(self.n_f_split[i]))
                else:
                    deepest_transit.append(900)
            tran_selec = numpy.argsort(deepest_transit)[0:ntransit]    
        
        lc = []
        t = []
        for i in tran_selec:
            lc = numpy.append(lc, self.n_f_split[i])
            t = numpy.append(t, (self.t_split[i]+self.porb*24*i)/24+self.x0)
        
        phase = (t % self.porb)/ self.porb
        jj = numpy.argsort(phase)
        ff = phase[jj]

        self.smoothed_LC = scipy.ndimage.filters.uniform_filter(lc[jj], size = 10) # equivalente ao smooth do idl com edge_truncade

        x = phase[jj]
        y = 1 - self.smoothed_LC
        yh = 0.002

        kk = numpy.where(y >= yh)

        x1 = min(x[kk])
        x2 = max(x[kk])
        fa0 = (x1 + x2)/ 2 # valor central dos transitos em fase

        self.time_phased = (ff - fa0)*self.porb*24

        bb = numpy.where((self.time_phased >= min(self.ts_model)) & (self.time_phased <= max(self.ts_model)))
        
        return self.time_phased[bb], self.smoothed_LC[bb]

    def gettime_phased(self):
        return self.time_phased
        
    #--------------------------------------------------#

    def select_transit_smooth(self, selection):
    
        '''
        Funcao para selecionar o transito desejado
        selection :: número do transito selecionado
        '''   

        i = selection
        
        lc = []
        t = []
        
        # for i in tran_selec:
        lc = numpy.append(lc, self.n_f_split[i])
        t = numpy.append(t, (self.t_split[i]+self.porb*24*i)/24+self.x0)
        
        phase = (t % self.porb)/ self.porb
        jj = numpy.argsort(phase)
        ff = phase[jj]

        self.smoothed_LC = scipy.ndimage.filters.uniform_filter(lc[jj], size = 10) # equivalente ao smooth do idl com edge_truncade

        x = phase[jj]
        y = 1 - self.smoothed_LC
        yh = 0.002

        kk = numpy.where(y >= yh)

        x1 = min(x[kk])
        x2 = max(x[kk])
        fa0 = (x1 + x2)/ 2 # valor central dos transitos em fase

        self.time_phased = (ff - fa0)*self.porb*24

        bb = numpy.where((self.time_phased >= min(self.ts_model)) & (self.time_phased <= max(self.ts_model)))
        
        return self.time_phased[bb], self.smoothed_LC[bb]